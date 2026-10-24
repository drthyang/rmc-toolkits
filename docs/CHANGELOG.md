# Development Log

Chronological record of notable changes, newest first. For current architecture and conventions see
[AGENTS.md](../AGENTS.md); for forward plans see [ROADMAP.md](ROADMAP.md).

## Unreleased

### Bond Geometry: the "Angles" tile is an angle count, and says so

- The KPI tile labelled "Angles" (shown uppercase, "ANGLES 15.0 per Ta") is a count, not an angle:
  `angleCount / apexCount`, the mean number of A–B–C angles at one B atom over every B atom of the
  configuration. It is now labelled **Angle count**; the hero card's ? help (now *About the plot
  and the angle count*) opens with its definition (all B atoms in the denominator, unbonded ones as
  0; one angle per pair of bonds, n(n−1)/2 for A = C; inclusive windows; every periodic image, no
  minimum-image cut-off; the one `.rmc6f` loaded, no averaging over configurations), and the tile's
  hover starts "Mean number of angles per B atom (a count, not an angle): N angles over M B atoms". The Model information card's lattice angles, which sat
  above it on the same page under the same word, read **Cell angles** (hover: α · β · γ).
- Display text only: the engines, the `/api/triplets` and worker payload keys (`angleCount`,
  `apexCount`), the CLI and the Copilot's `compute_bond_angles` data are unchanged. An independent
  brute-force count (scipy `cKDTree`) matches the engine's totals exactly on the Demo
  ([bond-geometry.md](algorithms/bond-geometry.md#step-3--the-kpi-rail)).

### Dashboard: the G(r) error ribbon, F(Q)'s σ propagated through RMCProfile's transform (protocol 2026-10-24.4)

- New module `rmc_toolkits/ft_uncertainty.py` (exported; the browser port is
  `src/ftUncertainty.js` with its run glue `src/ftUncertaintyRun.js`). It rebuilds the transform
  RMCProfile writes as `<stem>_FT_XFQ<n>.csv`: `X_ray_exp_renorm` = Σ ΔQ·Q·(F + c)·sin(Qr) over the
  run's F(Q) grid, a rectangle rule with no 2/π and no ρ₀ (Keen's D(r) × π/2), no Lorch, c the run's
  fitted offset. The run's `F(Q)_Expt` is the `.dat`'s input file renormalised by the fitted scale
  s and, with `CONVOLVE`, box-convolved at the PDFpartials' last r, minus c. Two reproduction tests
  establish the chain on the run's own files before anything is drawn, each to 10⁻⁶ of the
  function's maximum. The demo's G(r) reproduces to 1.5×10⁻¹⁴, and two local runs' F(Q) from their
  input files to 5×10⁻¹⁵.
- σ comes only from the input file's own third column, clean by `usable_sigma`'s rule, times s. It
  is propagated exactly: σ(r) and the full covariance ½[P(r − r′) − P(r + r′)], with the correlation
  between r points and the number of independent points. The box convolution is left out of the
  propagation and checked exactly to change σ(r) by under 1 % (0.19 % at r ≤ 20 Å on the demo's
  grids).
- The Dashboard's xPDF card draws the ribbon (InteractivePlot's opt-in `ribbon`: a band under the
  experiment, in the legend, ± σ in the readout). When the G(r) file has its own σ column, that one
  is drawn and the propagated one is a toggle. The caption says the r points are correlated, and
  its ? help gives the transform, both reproduction residuals and the tolerance. With no σ it says
  why: on the demo, "F(Q) carries no uncertainties: scale_ft_rmc.fq (named in GTS_250K.dat) is not
  in the run folder". A failure is also a Problems-section warning, and Live Data re-reads the
  ribbon with its files.
- Flask's `/api/plot/data` sends the same payload as `uncertainty` on an xPDF plot (cached on the
  transform and its companions; the input read only inside the data roots). A column naming an
  uncertainty is the experimental series' `sigma`. One golden holds both runtimes to 10⁻⁹:
  `ft_uncertainty_fixture.json` (`tests/generate_ft_uncertainty_fixture.py`). Known answers cover a
  brute-force Jacobian, Monte Carlo, white noise, the Nyquist grid, MATERIA's pystog golden adapted
  to RMCProfile's convention, and StoG (data-gated).
- `compare_fit` calls an X-ray G(r)'s residual feature above noise only against σ(r) (≥ 3 σ). Its
  `noise` block gives the σ's source, the residual's rms in σ, the reproduction residual and
  tolerance, and the correlation. A G(r) without σ(r) says so and calls nothing above noise.
  Model-visible: its description (+56 tokens on the worst turns; 48k headroom 5,200 dock and
  7,999 arm B) and the Demo's G(r) result (new snapshot `compare_fit.xpdf.json`). Arms A and A+ on
  dev are unchanged. Approved by the manager on the user's delegation, 2026-10-24
  (`validation/findings/decisions/2026-10-24-ft-error-propagation.json`).

### Copilot: `polyhedral_distortion`, a polyhedron's distortion indices (protocol 2026-10-24.3)

- A new engine, `rmc_toolkits/distortion.py` (`polyhedral_distortion`, exported from
  `rmc_toolkits`; source of truth) with its straight port `workers/distortion.js`: a centre
  element's coordination polyhedra — its ligands within an inclusive window, from the bond-angle
  engine's own neighbour search, the polyhedra of the most common count measured and the rest
  counted — by Baur's bond-length distortion index DI = (1/n) Σ |lᵢ − l̄| / l̄ (*Acta Cryst.* B30,
  1974) and, for octahedra and tetrahedra, Robinson's quadratic elongation ⟨λ⟩ = (1/n) Σ (lᵢ/l₀)²
  (l₀ of the regular polyhedron of the same volume) and bond-angle variance σ²θ over the cis
  angles (12 of an octahedron, the trans pairs being the matching with the largest angle sum;
  6 of a tetrahedron; θ₀ = 90° or arccos(−1/3)) (*Science* 172, 1971). Each is summarised (mean,
  sd, 10th, 50th and 90th percentile) over the configuration's polyhedra and over the average
  structure's (the site means of the average-structure basis), so a static distortion separates
  from thermal spread.
- The PCA worker's `distortion` request and `GET /api/distortion` (`rMax` ≤ 8 Å; a window holding
  more than 10⁶ bonds is refused before any is stored); the request contract gains the route, its
  limits and 17 cases. Known answers: regular polyhedra 0, 1 and 0; the 4 × 1.93 + 2 × 2.25 Å
  Jahn–Teller octahedron DI = 0.069831 exactly; trigonal and flattened polyhedra and an off-centred
  atom in closed form; a thermal-only network whose mean DI is √(2/π)·√(5/6)·s/l̄ for a bond spread
  s, and whose average octahedron is regular. A parity golden (`distortion_fixture.json`) holds the
  port to 1e-9, refusals included.
- The AI Copilot's `polyhedral_distortion(centre, ligand, r_min, r_max, window_confirmed?)` renders
  the configuration's and the average structure's indices. It refuses a window that cuts the
  coordination shell by `compute_bond_angles`' rule (the Workbench asks first), unless the user chose
  it. Its description says only what it computes, which is what the validation's arm C reads; the
  method advice (where the window comes from, octahedra or tetrahedra isolated or in an extended
  network, the static/local reading) goes to the dock and arm B. On the Demo the TaSe6 octahedra
  carry a static 3 + 3 distortion (average structure DI 0.0415, configuration 0.0469), while the GaSe4
  tetrahedra are regular on average (0.0010 against 0.0278).
- A fact sync of the interviewed skills, the tool being on main: coordination-shells step 7 names
  `polyhedral_distortion` and what it returns (it said no Workbench tool computes a distortion
  index), local-order step 4 says what it returns (it said "where the Workbench offers it"), and the
  skill lint, which now finds the tool among those the Copilot offers, drops its allow-list entry
  ("awaits the polyhedral_distortion PR").
- What the model sees: the tool's definition (+1,745 characters in the dock, +1,741 in arm B, +1,246
  in arm C), the first request's tool lists, a new tool-results sample, and read_skill's
  coordination-shells text (+160 characters). The 32k budget pins move by +812 (dock) and +808
  (arm B) tokens; at the live runs' 48k the dock keeps 4,873 tokens and arm B 7,607, with the
  interviewed skills. Arms A and A+ do not read the tool, and their dev records are unchanged.
  Approved by the manager on the user's delegation, 2026-10-24 (the user's interview answer 7 for
  coordination-shells; `validation/findings/decisions/2026-10-24-polyhedral-distortion.json`).

### Tests: the χ² clamp compared to 1e-12 relative (oldest-deps CI)

- `test_log_data_carries_every_chi_term_and_their_total` compared the zero curvature term's clamp,
  `np.log(1e-12)`, with `math.log(1e-12)` bit for bit. Under the oldest pinned numpy its SIMD log
  can differ in the last ulp, depending on the runner's CPU, which turned main's Python 3.9 job red.
  The test now checks the same things exactly — one value on every row, the right row count — and
  the clamp's value to 1e-12 relative. The code is unchanged.

### Copilot: the interview answers of four skills as checks — coordination shells, local symmetry, displacement anisotropy, local order (protocol 2026-10-24.2)

- **coordination-shells.** Every coordination number or bond angle in an answer carries its window
  in Å: an answer that gives one without it — in its sentence, or stated once for the pair or the
  whole answer — gets a `window` finding under it (the user's rule, 2026-10-24). The coordination
  check reads "4 F per Cu" as Cu's count of F, and a clause that states the whole shell's window as
  a claim about the shell (before, any window in Å made the clause partial). `find_gr_peaks` gives
  the average structure's first shell for the pair, `first_shell.average_structure_d_A` (its
  shortest and longest bond between the site means), and warns when a bond lies outside the window
  by more than the first peak's FWHM: no Demo or dev pair does (the Demo's Ta–Se 2.52–2.76 Å in
  2.21–3.23 Å); a 4 + 1 pyramid whose shell closes before its fifth bond does. The evidence gate
  backs distorted octahedra or polyhedra with `polyhedral_distortion` where the run offers it, else
  `analyze_clusters` (EVIDENCE_RULES says so), an `anyOf` row now asks for the first analysis the
  turn offers, and `analyze_clusters`' method keeps it to finite clusters. Approved by the user,
  2026-10-24 (`validation/findings/decisions/2026-10-24-coordination-shells-interview.json`).
- **local-symmetry.** A group–subgroup relation the answer claims that does not hold gets a
  `subgroup` finding: read against the run's own groups where it knows both (the ladder's rungs are
  nested; the index is |G|/|H| in the cell; t when the translations per cell stay, k when the point
  group stays, else a general subgroup, never maximal — the Demo's P-42₁m in F-43m is index 12,
  neither t nor k), and by the point groups alone for a published phase the run does not hold
  ("F-43m → P4/mmm" cannot be: 4/mmm has an inversion centre -43m lacks). A cell doubling the box
  cannot hold — odd along an axis it doubles — gets a `cell-doubling` finding, as does a doubling
  along axes no ordered zone-boundary point of `analyze_modulations` doubles. The symmetry-averaged
  CIF calls itself an "average-structure model of an RMC configuration", never a refined structure.
  Approved by the user, 2026-10-24
  (`validation/findings/decisions/2026-10-24-local-symmetry-interview.json`).
- **displacement-anisotropy.** `analyze_displacement_directions` labels its peak in the user's
  three tiers: a **finding** at ≥ 3σ (unrounded) and ≥ 1.2× isotropic that holds at the
  neighbouring resolution (ν − 1; ν + 1 at ν = 1), which the tool reads for a candidate — the peak
  test passes there and the two peaks lie on one axis within the two cells' radii; a call at
  Auto ± 1, as the skill asks, reads the Auto resolution as the other; a **hint** at
  ≥ 3σ under 1.2×, at 2.5–3σ or not held; **nothing** below 2.5σ. The second reading costs 0.7 ms
  on the Demo and about 240 tokens per finding. On the Demo, Ta site 17 and Se site 34 are
  findings and Ta site 25 (4.47σ at ν = 2, 2.33σ at ν = 1) is now a hint; every dev split site is
  a finding. The direction check and arm A+ take a finding only (no dev verdict moves). An answer
  that calls a split or flat-topped site disorder before `analyze_modulations` tested its element
  gets a `disorder` finding (ordered antiphase displacements give the same split cloud), and the
  evidence gate reads a flat-topped site called disordered as it reads a split site. Approved by
  the user, 2026-10-24
  (`validation/findings/decisions/2026-10-24-displacement-anisotropy-interview.json`).
- **local-order.** "No order" now also needs `polyhedral_distortion` on each centre of octahedra or
  tetrahedra no finite cluster holds (the Demo's Ga, a perovskite's B site), where the run offers
  it; its diffuse verdict is `analyze_modulations`' own. The diffuse verdict is read as the user
  reads it: only antiphase (or general-q) order past its bar is order; in-phase correlation is
  "neighbouring cells moving together (correlated motion), not order", as `analyze_modulations` now
  says it, and never supports "no correlation". Three checks on the answer: `no-order` ("no order"
  when order was found), `correlated-motion` (order between cells when only in-phase correlation, or
  none, was found) and `no-correlation` (uncorrelated displacements when a correlation passed, in
  phase included). Arm A+ claims an antiphase verdict as a modulation with no q, the verdict's one
  kind for order between cells; no dev case passes the diffuse test, so no dev record moves.
  Approved by the user, 2026-10-24
  (`validation/findings/decisions/2026-10-24-local-order-interview.json`).

### Copilot: the fit-assessment interview's rules in code — convergence on every χ² term, the 5 % constraint, two answer checks (protocol 2026-10-24.1)

- **Every χ² term of the log and their total.** Both runtimes read every χ² column after `Time
  moves_acc moves_gen` (`readChi` / `read_chi_log`, `CHI_TERMS_START`), each row's total and the rows
  per log, and hand them to the Copilot and the watchdog as `plotData.chiTerms` (`/api/plot/data` sends
  them in Flask mode); the per-term histories concatenate only when every log names the same terms. The
  Dashboard's chart still draws the last column.
- **Converged only when no term and not the total still falls beyond its noise** over the last tenth
  of the run (the user's rule, 2026-10-24). A series falls when the least-squares line through its last
  tenth (at least 10 rows, never back past the last restart) drops by more than 3 × the scatter of its
  points about that line — the app's 3σ bar: a flat series of independent noise passes for falling
  0.8 % of the time over 10 rows and never over the Demo's 227 (`src/llm/watchdog/heuristics.js`:
  `seriesTrend`, `judgeConvergence`; the old fixed 0.02 in ln χ² is gone). `assess_convergence` reports
  every term and the total (change, noise, trend), the context's `convergence.last_10pct` replaces
  `change_over_last_10pct` and `recent_slope_per_step`, the system prompt says what `watchdog_status`
  judges, and the watchdog badge moved to the same rule (AGENTS.md next step 6). The Demo stays
  `improving` (F(Q)_1, X_ray_(R)1 and the total fall; its curvature constraints are flat), as does
  dev-fit-unconverged.
- **Constraint-limited above 5 %.** `find_gr_peaks` calls a peak constraint-limited when g(r) just
  above the closest-approach constraint is above 5 % of the peak (`CONSTRAINT_LIMIT_FRACTION`; the
  shell separation keeps 2 %), and advises interpreting nothing on its low-r side and rerunning with
  the constraint loosened (method: arm C reads the reading alone). No reading moves on the Demo
  (Ga–Se 8.7 % at 2.2 Å) or the dev cases (no constraints).
- **Two checks on the answer:** `long-wave` flags order claimed at a box-scale wavevector (every
  |q_i| ≤ the modulation engine's 0.25, or a long wave in words): a long wave says nothing about order
  at its own q; `single-configuration` flags a distribution's numbers (bond angles, coordination
  counts, displacement directions, cluster shapes) quoted without saying they come from one
  configuration. On the answer corpus they raise one finding, on the README GIF's GaSe₄ answer
  (judged right), and no false alarm. Scripted replays on the Demo pin the user's answers 1, 2, 5 and
  6 wherever a check can catch the wrong answer.
- The skill lint binds the interviewed rmc-fit-assessment skill's four lines this code builds — the
  every-term convergence rule (with a probe), the 5 % gate to `CONSTRAINT_LIMIT_FRACTION`, and the
  two firm rules to the `long-wave` and `single-configuration` checks (with probes) — so their
  allow-list entries ("awaits the fit-assessment code PR") leave.
- What the model reads: the system prompt +226 characters; `assess_convergence`'s description 274 →
  575 characters and its Demo result 435 → 1,111; the Demo context +486 characters, with the context
  budget raised from 6,000 to 6,500 characters so its second g(r) peaks stay. Arms A, A+ and D on dev
  are unchanged. Approved by the user, 2026-10-24
  (`validation/findings/decisions/2026-10-24-fit-assessment-interview.json`).

### Skills: the facts synced with the tools (the user's decision of 2026-10-19)

- The batch-1 skill drafts now state the numbers and rules that the tools apply. Nothing else in
  them changed: no method, step order or judgment call.
  - The Demo's Ta–Se first peak is 2.21–2.65 Å, holding 3.25 Se per Ta. The drafts said
    2.33–2.71 Å and 3.9, from before shells.js moved the first peak to the inter-peak minimum.
  - A split site is the app's rule: a resolved PC1 whose excess kurtosis is ≤ −0.5 and at least 3
    sampling errors below zero. A negative Mardia kurtosis alone is not one.
  - A preferred direction needs the peak test at 3σ or more on the unrounded σ.
  - The low-index cone is "within 15°", 15° included. That is under half the 35.3° between <110>
    and <111>, not half of it.
  - A family is named only when the peak cell's radius is 15° or less (plan item 2.5).
  - compute_bond_angles' cut rule is the CN margin (max(0.5, 15 % of the shell's CN)), as #113
    made it. The counts are its fallback when there is no shell to compare.
  - analyze_clusters runs for elements whose own first shell holds 0.5–4.5 neighbours of their
    own kind. It reads the shape against a shuffle null at each cluster's own size, and the size
    against independent atoms.
  - The average structure's group is read from 1.25 × the noise floor.
  - The order rule is ratio ≥ 3 and ≥ 3σ, against independent cells (F1).
  - The cluster tool's z is no tail probability.
  - The Displacement Directions page opens at Auto, the engine's resolution, with its numbers from
    the raw histogram (the skill said ν = 10).
  - A wrong reason in rmc-fit-assessment step 7 (long waves "left out of the background") is
    removed.
- The skill lint's allow-list holds no FACT entry any more, #113's three "awaits the skill
  fact-sync PR" entries included. One JUDGMENT stays: how to word a rung just above the noise
  floor. The bindings follow the new text:
  - the CN margin is bound to `CN_CHANGE_ABS` and `CN_CHANGE_REL`, with #113's
    compute_bond_angles asserts;
  - the counts fallback to `MODAL_COORDINATION_SHARE` and `UNBONDED_APEX_FRACTION`;
  - the split rule to `SPLIT_SITE_KURTOSIS` and `SPLIT_SITE_SE`;
  - the cell radius to `directionResolved`.
- `docs/SKILL_INTERVIEWS.md` lists the eleven choices the sync found and left to the user. The
  read_skill and over-claim record snapshots are re-recorded. The 32k budget pins move by 853 tokens
  (dock) and 854 (arm B). At 48k, 6,886 and 9,615 tokens are left to spare.
- **rmc-fit-assessment interviewed (2026-10-24).** The skill is now the user's own answers. Its
  status reads "Interviewed 2026-10-24", and `docs/SKILL_INTERVIEWS.md` §1 records each answer and
  where it landed. The rules it now states:
  - convergence judged on every fit term and the total;
  - a "still improving" caveat only for first-shell quantities;
  - no fixed Rw cutoff;
  - a peak constraint-limited above 5 %;
  - at least 3 independent runs and 10 accepted moves per atom past convergence before a
    distribution;
  - the user's list of what is never a finding, and the report's always and nevers;
  - two new firm rules;
  - the user's wording for long waves.

  The code half (every-term convergence, the 5 % constant, the two answer checks, the replays) is
  its own PR. Until it lands, four skill-lint entries are marked "awaits the fit-assessment code
  PR". `copilotSkills.test.js` tells an interviewed skill from a draft.
- **coordination-shells interviewed (2026-10-24).** The skill is now the user's own answers, and
  §2 records where each landed. The rules it now states:
  - the window cross-checked against the average structure's bond lengths;
  - a split shell's long − short ± its sampling error always quoted;
  - a shell not separated from the next cut at the minimum, with the 0.2 Å range beside it;
  - the cut rule kept, with no exemption for disordered or glassy systems (the user's window
    instead);
  - angles read from the sin θ-corrected histogram;
  - a distortion needs the cluster effect at z ≥ 3, the same sign in ≥ 3 independent runs, and
    angle spreads;
  - extended networks get angle spreads plus a polyhedral distortion index;
  - a new firm rule: every CN or angle carries its window. It waits for a later coordination code
    PR, as a skill-lint entry.
- **local-symmetry interviewed (2026-10-24).** The skill is now the user's own answers, and §3
  records where each landed. The rules it now states:
  - the 1.25 × margin, with a rung below it "at the limit of the sampling", never a lowering;
  - "the" group is the lowest rung at or above the margin;
  - a lowering is described by its lost operations plus the group–subgroup index and its t or k
    type;
  - against a published cell-doubled phase, only what the box can represent is compared;
  - CIF export rules;
  - groups named as "consistent with ⟨group⟩ at ⟨tolerance⟩ Å";
  - two new firm rules, which wait for a later local-symmetry code PR as skill-lint entries.

  The lint's last JUDGMENT entry is gone: the margin is bound to `NOISE_FLOOR_MARGIN`.
- **displacement-anisotropy interviewed (2026-10-24).** The skill is now the user's own answers,
  and §4 records where each landed. The rules it now states:
  - "the most disordered site": split sites first, then U_iso and anisotropy;
  - a negative-kurtosis site that misses the split rule is "flat-topped, not split";
  - an element pooled only for an element-wide summary, with the frames caveat;
  - three tiers for a direction: a finding at ≥ 3σ and ≥ 1.2× isotropic, a hint at ≥ 3σ under
    1.2× or at 2.5–3σ, nothing below 2.5σ;
  - a finding holds at the Auto resolution and at a neighbouring ν;
  - every reported direction carries ± its cell's radius, and the 15° cone stays;
  - a cluster distortion needs z ≥ 3 with the same sign in at least 3 independent runs;
  - a new firm rule: a split or flat-topped site is called disorder only after
    `analyze_modulations` ran for its element.

  The tiers, the two resolutions and the firm rule wait for a later displacement-anisotropy code
  PR, as five skill-lint entries. The skill is now the largest the Copilot reads (8,896 of
  read_skill's 9,000 characters). The 32k budget pins move by 868 tokens, and at 48k 5,685 (dock)
  and 8,415 (arm B) tokens are left to spare.
- **local-order interviewed (2026-10-24): all five skills are now the user's own.** §5 records
  where each answer landed. The rules it now states:
  - the list adds `analyze_modulations`' diffuse verdict and, for an extended network,
    `polyhedral_distortion`;
  - the diffuse verdict, past its calibrated bar, counts as order between cells;
  - directions for every symmetry-inequivalent site, not only the outliers;
  - a partial under 5 % of the signal is "not constrained by the data", and "no order" leaves it
    out;
  - never a bare "no order": a test that could not run is named with its reason;
  - a new firm rule: "no order" is refused until every applicable analysis has run.

  Three skill-lint entries wait for their code: the extended evidence gate and the diffuse
  verdict's reading (a later local-order answer-check code PR), and `polyhedral_distortion` (its
  own PR). No snapshot or budget pin moves.

### Validation: the held-out marginal pair re-set (protocol 2026-10-19.20)

- The marginal situation was re-set after F1's corrected null raised the test's power. A fresh author
  who had seen no dev value re-chose the two held-out marginal matches' amplitude from the updated
  `judgment-spec.md`, with new seeds by a rule fixed beforehand; the values are sealed
  (`TEST_SET_SHA256` changes). The lead checked the author's blindness from its tool log. The holdout
  scan no longer reads the digits that end an identifier (`SHA256`) as a number. Arms A and A+ on the
  test set are re-recorded with `--final` from a clean tree. Approved by the manager on the user's
  delegation, 2026-10-19; the user may reverse it before the freeze.

### Validation: the marginal situation re-set for the independent-cell null (protocol 2026-10-19.19)

- **Why.** The modulation engines' calibrated null (F1, protocol 2026-10-19.16) gives an isotropic
  site three components where the earlier null took one. The dev weak-order pair's amplitude had
  been set from that one-component bound to an expected 3σ, so both of its seeds read far above the
  bar, and "a signal at the edge of significance" was no longer at the edge.
- **The spec** (`validation/judgment-spec.md`, the held-out author's only description of the knob)
  gives the expected significance under the new null:
  - the B power at q is σ²/(3σ² + δ²) · χ'²₃(nδ²/σ²);
  - the null is Σ λ_k χ²₁ with λ = n/(n − 1) · (σ², σ², σ² + δ²)/(3σ² + δ²);
  - the Šidák family is every zone-boundary point of every element;
  - the expected significance is the median realised σ, the σ of the median power.

  The spec also says how to evaluate it and how to check it with the engine on fresh
  configurations. It holds no dev value, and `judgment.test.js` still scans it.
- **The dev pair.** Its amplitude is re-set by that formula to an expected 3σ, and the engine over
  fresh configurations agrees (about half clear the bar). Its seeds come from the pair's own rule,
  fixed before any σ was read at the new amplitude; every other parameter is kept. The values are
  in `validation/findings/decisions/2026-10-19-marginal-reset-dev.json`, which the held-out author
  does not read. `judgment.test.js` checks the formula with the engine's own functions, the pair's
  expected 3σ and its seed rule.
- **Arm C.** `analyze_modulations`' long-wave and in-phase notes lose their causal parentheticals
  ("(common in RMC)", "(as RMC's long waves and acoustic-like thermal motion give)"): what a pattern
  is evidence of is inference, which arm C drops while keeping definitions and limits. Both notes
  are built with `methodCheck`. Arm B and the dock read them byte for byte as before, and every
  model-visible snapshot passes unchanged. New arm C samples:
  `tool-results/analyze_modulations.method-off.json` and `analyze_modulations.element.method-off.json`.
- **Arms A, A+ and D on dev.** Only the two re-set cases' six records move. One A+ record still
  claims the weak order; the other, below the bar, no longer does. They are re-recorded in
  `validation/baseline/arm-a-dev.json`; every other record and every arm A and D score is unchanged.
  The held-out marginal matches are re-set from the spec by their author, separately.
- Approved by the manager on the user's delegation, 2026-10-19; the user may reverse it before the
  freeze.

### Copilot: one cut rule — compute_bond_angles refuses on the Workbench's CN margin (protocol 2026-10-19.18)

- Whether a bond window cuts its coordination shell is decided once, by the shell module's
  `windowCut` (`src/llm/context/shells.js`), for the Workbench's question before an analysis runs and
  for `compute_bond_angles`' refusal: a window that meets its first shell is cut when it holds fewer
  partner atoms than the whole shell by at least max(0.5, 15 % of its CN), whatever its coordination
  counts. Until now the tool refused on the counts (more than half the central atoms unbonded, or
  fewer than 60 % with the most common count) while the Workbench asked on the CN margin; the counts
  now decide only where the run reads no shell to compare the window with: no partials file, no
  closed shell for the pair, no density for a CN, or a window apart from the first shell (the module
  reads first shells only). Approved by the user, 2026-10-19
  (`validation/findings/decisions/2026-10-19-shell-cut-rule.json`).
- What changes: a whole shell is never refused (dev-offcentre's Nb–Nb, which runs into the next
  shell, 52 % of its Nb at the modal count), nor a window within the margin (the top 0.2 Å early);
  a sub-shell first peak is refused unless the user chose it (the Demo's Ta–Se, 3.25 of 6;
  dev-jahn-teller's four short Cu–F). On the 149 natural windows of the decision record, 15 change
  verdict (11 no longer refused, 4 refused now), and the Workbench asks on every window the tool
  refuses (`validation/run.js shell-cuts` shows both beside the old counts rule).
- What the model reads: the refusal names the CN shortfall by path ("it holds 3.75 Se per Ta
  (shell_ab.cn_in_window) of the shell's 6 (shell_ab.cn, in shell_ab.window_A)…"); every result adds
  the shell reading (`shell_ab`, and `shell_bc` for a B–C window of its own) and the counts' `modal`,
  `modal_share` and `unbonded_share` as information, with no note for a count spread; the
  coordination answer check counts a window as the whole shell when its CN is the shell's within
  the margin. The tool's description and schema are unchanged; five `compute_bond_angles` tool
  results and four Workbench tool messages grow by about 200–250 characters.
- `buildCopilotTools` takes the run's `plotFiles` and `runSettings`, so the tool reads the shell
  catalogue the Workbench reads; `window()` loses its `counts` option and the llm module exports
  `windowCut` instead of `countsCut` and its two thresholds, which stay in the shell module as the
  fallback.
- The coordination-shells skill still states the counts rule: the skill lint binds its 60 % and
  "most" to the CN margin and holds the three lines on its allow-list, with both sides, until the
  skill fact-sync PR. Arms A, A+ and D on dev are unchanged (A+ reads no bond-angle window).

### Displacement Directions: the page and the tool at the same settings, and the direction's resolution (plan item 2.5; protocol 2026-10-19.17)

- **The page opens at Auto**, the engine's recommended resolution (the finest ν whose cells still
  average at least 12 displacements), which is the tiling `analyze_displacement_directions` reads.
  Until now it opened at ν = 10. The Auto choice is labelled with the ν the engine chose for the site
  on screen (`Auto · ν=2 (42)`), also while another resolution is picked. A picked resolution
  overrides Auto and stays, like every option on the page: across a site change, a tab switch and
  another run, until the user picks again or a Copilot step's Show in sets it.
- **Every readout comes from the raw histogram.** Both engines pick the peak cell on the smoothed
  map, so with smoothing on, the peak's cell, direction and σ moved with the Smoothing slider (on a
  dev case at ν = 10, two passes turned a −0.28σ peak into −10.4σ). The page now asks for the
  readouts at smoothing 0 (the tool's request) and draws the same request smoothed: the smoothing
  only redraws the map, its colorbar and the hover's enhancement. A smoothing change asks for the
  map only. The map test, the ± asymmetry and the Bingham anisotropy never depended on it.
- **The map opens unsmoothed** (smoothing 0, was 2). Auto's cells already average at least 12
  displacements, and on its coarse tilings two passes left about 30 % of a cell's excess in place.
  This is UI only: no readout, and no score of the validation's arms, depends on it.
- **The direction's resolution is reported.** The peak direction is the peak cell's centre, so it
  names a direction only to within the cell. Both engines (`orientation.py`, `workers/orientation.js`)
  now return `peakCellRadiusDeg`, the peak cell's circumscribed radius: the largest angle from its
  centre to its corners, so every direction counted in the cell lies within it — 37.4° at ν = 1,
  19.3–20.9° at ν = 2, under 15° from ν = 3. The page prints the peak "± R°, the cell radius at
  ν = n"; `analyze_displacement_directions` gives `peak.radius_deg` and the same wording, and above
  15° (`LOW_INDEX_CONE_DEG`) says the direction is not resolved finer than R° and names no <100>,
  <110> or <111> family. The orientation golden pins the field in both runtimes.
- One module holds the page's defaults and both requests (`src/directionSettings.js`), and
  `src/useDirectionHistogram.js` asks for both; `OrientationView` draws `result` and prints its
  readouts from `readout` (`DirectionReadouts`).
- **Validation:** arm A, and A+ through the same reading, reads a split site's direction map at
  the page's own defaults (`ARM_A_RULES.DIRECTION_SETTINGS` is `DIRECTION_DEFAULTS`); the reading
  records the ν it read and the peak cell's radius. Arm A takes a significant peak whatever its
  radius, as a user reading the page; arm A+ applies the tool's bound and takes PC1 where the peak
  cell is wider than 15° (`SWEEP_RULES.directionResolved`). `PROTOCOL_VERSION` 2026-10-19.17. On the
  dev set nine A and A+ records move, every one through that reading. One judgment case's arm A
  score falls: a significant peak on Auto's 12-cell tiling names a direction 31.7° off its cube axis
  (beyond the 15° tolerance) where PC1 named it before. Its A+ score is unchanged, since A+ keeps
  PC1 there. The moved records are re-recorded in `validation/baseline/arm-a-dev.json`, the others
  as recorded. The test set was not run.
- Known answers: on every site of the dev thermal case the page's default readout equals the
  tool's result field for field, and arm A reads exactly that; its flattened O sites clear 3σ on
  the map test there and stayed below it on the old default. A planted [0 0 1] split site of 216
  copies reads ν = 1, ± 37.4° and the unresolved wording; at 4000 copies, ν = 5 and resolved.
- What the model sees changes in `analyze_displacement_directions` only: its description names the
  cell radius (+21 characters a request) and its result gives `peak.radius_deg`, a longer note and
  the new wording. The model-visible snapshots, the record snapshots and the context-budget pins are
  updated for it. The displacement-anisotropy skill still says the page defaults to ν = 10 (two
  sentences), for the skill fact-sync PR.
- Decision record: `validation/findings/decisions/2026-10-19-directions-auto-resolution.json`.
  Approved by the manager on the user's delegation, 2026-10-19.

### Modulations: the independent-cell null (F1) and diffuse order between cells (plan item 4.4; protocol 2026-10-19.16)

- **F1.** The modulation engines (`modulations.py`, `workers/modulations.js`) test an element's
  power against its own cells placed at random: every correlation inside a cell is kept (the
  ligands of a Jahn–Teller octahedron, the vertices of a cluster), every correlation between cells
  destroyed. The null's cross-spectral matrix is white, N/(N − 1) Σ_c y(c) y(c)ᵀ at every q ≠ 0, so
  the power is a weighted sum of χ²₁ (zone boundary) or χ²₂/2 (general q) terms with its
  eigenvalues as weights (`independent_cell_null` / `independentCellNull`); the tail is a
  Lugannani–Rice saddlepoint, 0 to +8 % above the exact tail from 10⁻² to 10⁻⁹ (`null_tail` /
  `nullTail`). Each element reports `nullDof`. The null of independent sites it replaces found
  order in 8.8–31.8 % of fresh configurations of the dev Jahn–Teller, two-feature and octahedral generator
  kinds (1,000 each); the calibrated null in 0.1–0.6 %, each of its three tests at its 0.135 %. It
  also stops holding an independent site to one component: a weak order on one isotropic site is
  found in 99 % of configurations where it was found in 44 %, and both dev weak orders are now far above the
  3σ bar.
- **4.4.** Each element gains `diffuse` and `verdict.diffuseOrder`: the ranks of the power over the
  compared modes (no long wave, no zone-boundary point, none found ordered) regressed on the 13
  neighbour-cell lags cos(2π q·R), χ²_r under independent cells, Šidák over the elements; diffuse
  order needs 3σ and a neighbour-cell correlation |Γ| ≥ 0.1. Its Finding carries one Evidence
  record in the cluster engine's form, `clears` recomputable. Parseval fixes the total power, so
  the test reads how the power is spread. Planted short-range order (neighbour correlation ∓0.27)
  is found in every configuration, in antiphase or in phase; no configuration without order between cells
  reads as it beyond its bar (`tests/test_modulation_calibration.py`: a fast CI subset, and 1,000
  configurations per family with `RMC_TOOLKITS_FULL_CALIBRATION=1`). One strong mode never makes it. The
  Demo's three elements are correlated in phase beyond the long waves.
- `analyze_modulations` shows a `diffuse` row [lag, correlation, σ] per element (F2's rule on the
  correlation's 0.1 bar), says antiphase correlation as short-range order and in-phase correlation
  in one note as neighbouring cells moving together, and reads its σ "against independent cells".
  A nearest miss at σ = −∞ reads "< −8σ" (it printed "nullσ"). Its result budget is 4,500
  characters. What the model sees changes in this tool's result only: both tool-result snapshots
  (`analyze_modulations.json`, `analyze_modulations.element.json`).
- Arm A+ reads cellDoubling, generalOrder and longWave only (no feature kind for diffuse order).
  On the dev set one A+ record moves: a weak order now above the bar is claimed (re-recorded in
  `validation/baseline/arm-a-dev.json`); every arm A and D record and every other score is
  unchanged. The marginal judgment situation needs its amplitude re-set from the calibrated null,
  the user's call. Approved by the manager on the user's delegation, 2026-10-19
  (`validation/findings/decisions/2026-10-19-modulation-null.json`, `2026-10-19-diffuse-order.json`).

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

### Shells: a joined shell's split by number (plan item 4.2, protocol 2026-10-19.15)

- The shell catalogue (`src/llm/context/shells.js`) gives a shell joined across a deep trough its
  `split`: the parts cut at the minimum of the smoothed g between each two successive resolved
  peaks (the first is the first peak alone), the last ending at the shell's boundary; each part's
  CN over exactly its window (the parts sum to the shell's), its share and its mean distance and
  spread from r²g; `longMinusShort`, the last part's mean minus the first's, with
  `longMinusShortSe`, the configuration's sampling error (each part's spread over √ its bonds in
  the box, the bonds taken as independent; not the fit's uncertainty); and `direction`, the count
  the join decided by. A 2 + 2 + 2 octahedron keeps three parts and long − short is the longest
  minus the shortest: in the E⊗e Jahn–Teller problem that measures the amplitude whatever the mix
  of the Q₂ and Q₃ modes, and the shares show the mix. A shell whose peaks overlap without a deep
  trough has no split (its minimum cuts through both bond groups).
- Known answers (`shells.test.js`): 4 × 1.93 + 2 × 2.25 Å at σ 0.05 Å on 0.005–0.02 Å grids reads
  4 + 2, shares 2/3 and 1/3, long − short 0.32 Å within 0.002 Å; 2 + 2 + 2 three thirds and 0.37 Å;
  the reverse direction; the cut's bias (+0.001–0.0017 Å at a 0.07 Å spread, +0.006 Å at the edge
  of the join); and the sampling error's 95 % coverage, 0.953 over 300 synthetic configurations.
- `find_gr_peaks` returns it as `first_shell.split` (`parts` with window, mean, sd, CN in the
  join's direction and share; `long_minus_short_A` ± `long_minus_short_se_A`), declares it in its
  dataShape, names it in its description, and its joined-shell note quotes the parts, long − short
  and the shares by path instead of the peaks' separation, so the model quotes the split by
  reference. dev-jahn-teller's Cu–F reads 3.98 F per Cu at 1.929 Å + 2.01 at 2.251 Å, long − short
  0.321 ± 0.0043 Å (planted 0.32). What the model sees: the three tool-list snapshots (+123
  characters each; the budget pins +48 and +49 tokens); the Demo joins no shell.
- Arm A+ takes a Jahn–Teller claim's magnitude from the catalogue's split where it joins the
  claimed shell (else the g(r) maxima's separation arm A reads). On the dev set one record moves:
  dev-jahn-teller's A+ claim reads 0.321 Å where it read 0.31, with the same score (re-recorded in
  `validation/baseline/arm-a-dev.json`). Approved by the manager on the user's delegation,
  2026-10-19 (`validation/findings/decisions/2026-10-19-jt-split.json`).

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
- Above `MAX_BASIS_SITES` = 2000 sites, the space-group finder's cap (`MAX_SYMMETRY_SITES`, pinned
  equal in both test suites), nothing in the app can use the basis, so the route sends `basis: null`
  and `basisOmitted: {sites, cap}`; a one-site-per-atom box of 52 000 atoms would otherwise send
  15 MB on every request. The page turns that count into the browser's own "not analysed"
  description: the card, its missing ladder and CIF button, the Copilot's `symmetry` context and its
  two symmetry tools (offered, answering "Not analysed") are the same as for the same file picked in
  the browser. Two things the browser derives from the rows are not sent: the sampling noise floor
  (the context's `noise_floor_A` and the lowering tool's default tolerance; null in the browser too
  when every site has one copy, the usual case above the cap) and the pooled-site count in
  `analyze_displacement_directions`' note for an element.
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

### Copilot: result statuses on the tool contract; the event stream is the transcript (architecture candidate #11)

- `src/llm/copilot/tools.js` decides how every call ends, and the agent loop only forwards it. A
  call ends as one ToolResult — `{ resultStatus, args, data, checks, error, facts?, ms, cached, note,
  choice? }` — with `resultStatus` one of `RESULT_STATUSES`: `ok`, `refused` (the tool's own `fail`
  check, no error), `error`, `unknown_tool`, `skipped` (the loop's refusals, a closed set of
  reasons: the calls per round, one question at a time, the questions per request). The
  constructors (`executeToolCall`, `skippedCall`, `unknownToolCall`, `repeatedCall`, `withNote`,
  `questionResult`, `toolResult`) check the invariants and throw a TypeError on a result that
  breaks one; `severityOf` is the old `stepStatus`, `ran` the old `executed`; `toolMessage` and
  `compactToolMessage` moved there verbatim.
- `tool_end` carries the result: `resultStatus`, `note` (what the app added for the model; for a
  Workbench question, which sends no message, the text the model read in its place) and `choice`
  (a question the user answered: `{ source, index, confirmed }`, `confirmed` the clarifier choice's
  flag, which replaces reading "Keep …" off the label), beside the severity. No store: the dock and
  the validation harness read the same events, and the verification records the ToolResult
  itself (a `refused` step still counts as an analysis that ran).
- The dock: each step keeps the status, the note and the choice; the question card and the screen
  reader's line read a question's state from them in one wording (a skipped question is announced
  "Skipped: the Copilot uses its default, …", as the card says it, no longer "Answered: …"; a
  question that never reached the user "Not asked: …"); `ToolStep` shows the note sent to the
  model.
- The validation record: each step gains `result_status`, `note`, `choice`, `ref` and
  `analysis`; a step running when a turn stops is listed as stopped; `failed_calls` and
  `unknown_tool_calls` are read from the status (unchanged values), and `refused_calls` is added.
  The transcript prints a non-ok status and the note, and labels a Workbench question's card as not
  sent to the model.
- Protocol 2026-10-19.13: two stricter readings, approved by the manager on the user's
  delegation. "Asked before" in `unsafe_actions` counts only a question that reached the user (not
  a malformed or refused one, nor a window the Workbench settled with no card); number provenance
  credits only what the model was told (each tool message it read, its error and note included;
  never a Workbench question's card). Decision records `2026-10-19-asked-before.json` and
  `2026-10-19-provenance-told.json`.
- Unchanged: every byte, note and order the model reads (every model-visible snapshot, and a new
  corpus of raw tool messages written first, on the code before the change); `isAnalysis`, the
  evidence gate, the Outcome and the computed verdict; arms A and A+; the replays; the answer
  corpus; the scripted records (pinned as files before the change, equal after it with the added
  keys left out). Tests: `tools.test.js` (every status, the invariant table), `toolOutcomes.test.jsx`
  (one turn through the dock and the harness: the drift list), `toolMessages.test.js` (the raw
  messages), the harness's known-answer runs for both readings. Details: [the llm
  README](../web_app/frontend/src/llm/README.md), *How a call ends*.

### Shells: a split M₃/M₄ cluster shell does not join (F5, a known limit)

- A split first shell joins only when its summed CN is an ideal polyhedral count (4, 6, 8, 12), so a
  metal triangle or tetrahedron whose edges a distortion splits stays cut at its first trough. Accepting
  CN 2 and 3 was tried and left out: on the Demo and the dev cases it helps no shell and breaks one
  (dev-cell-wrap's V dimers read as an extended network; arm A+ there 1.00 → 0.83). The limit is
  documented with its evidence in `docs/algorithms/ai-assistant.md` and pinned by a known-answer test.
  Decided by the manager on the user's delegation, 2026-10-19
  (`validation/findings/decisions/2026-10-19-f5-shell-join.json`); the user may reopen it.

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

### Copilot: one turn policy decides what ends a turn (architecture candidate #10)

- `src/llm/copilot/turnPolicy.js` owns every counter that ends a turn — the rounds of tools, the
  turn's time, the rounds that added nothing, the one reminder and the warmer retry, the
  last-request note, the answer the evidence gate withdrew — and the order of the steps, behind one
  event-driven interface: `createTurnPolicy({ limits, clock, tools, judge }).next(event)`, the
  events `start`, `reply`, `round` and `no_tools`, the steps `ask` (why, tool choice, retry, notes,
  the verification's send-back), `run` and `end` (the `done`, and the text to say first). The agent
  loop reports each event and does what the step says; the policy judges each answer through the
  verification's `judge` and decides when (a send-back only with tools on, before the last request
  and with time left; the `stoppedBy` an answer is judged under; the withdrawn answer stands), the
  verification what. Turn time is read at most once per event, never in round 0, and the user's
  waits stay off it.
- One `done` shape on every path, built only by the policy: `{ type, text, rawText, outcome,
  rounds, stoppedBy, cut, reason, audit, references, gates }`. An answer gains `reason: null` and
  always has `rawText`; the "budget used" answer gains `rawText`, `references: null`, `gates: []` and
  `reason: null`, its Outcome kept to status, note and body (adding `computed` / `reasons` would
  change the harness's `record.final`); an empty reply gains `rawText: ''`, `references: null` and
  `gates: []`.
- `STOPPED_BY` is the one set of endings — answer, max_rounds, no_progress, empty, stopped, error —
  each with `cutShort` and the dock's note, and `stoppedByFailure` maps a thrown turn to 'stopped'
  (an AbortError) or 'error' (a stall included). The turn's `stopped` / `failed` events carry it,
  the dock closes a turn on it, `CopilotView` shows the row's note (its STOP_NOTES are gone, the
  sentences unchanged), the computed Outcome reads `cutShort`, and the validation harness records
  `record.stopped_by` from every terminal event.
- What left the agent loop: its rounds, start time, waits total, no-progress count, stall flag,
  reminder count, warmer flag, last-note flag and `why`, both turn-time tests, the stoppedBy
  derivation and the two inline `done` literals; the notes and fixed answers moved verbatim. The
  calls per round and the questions stay in the loop: they refuse a call and never end a turn.
- Unchanged: what the model sees (every model-visible snapshot), every scripted replay, the
  agentLoop and robustness assertions, the answer corpus, the harness's expectations and its
  `record.final` on every path (checked against the previous code on 23 scripted endings in both
  arms), the turn tests. No protocol bump: what a run sends is unchanged; `record.stopped_by` is a
  new, derived field. Tests: `turnPolicy.test.js` (rows of events → steps, the policy in the loop)
  and `src/__tests__/turnEndings.test.jsx` (every ending through the dock and the harness alike).
  Details: [the llm README](../web_app/frontend/src/llm/README.md), *The turn policy*.

### Modulations: the polarisation of an ordered mode (plan item 4.3), F2 and F3 (protocol 2026-10-19.12)

- `displacement_modulations` (`modulations.py`, `workers/modulations.js`) gives every zone-boundary
  point and strongest mode its `polarisation`: the leading eigenvector of the real part of the
  element's cross-spectral tensor at q, (1/m) Σ_s Re[U_s U_sᴴ] / (n_s tr C_s), whose trace is the
  power. It comes as a Cartesian unit `direction` and a direct-lattice `uvw` (largest component 1,
  the form of `list_sites`' `pc1_uvw`), both in the sense that makes the closest crystal axis acute,
  with its `share` λ₁/tr (1: one direction; ⅓: none preferred). The share is not a significance:
  one site at a zone-boundary point gives share 1 by construction. `element_spectrum` /
  `elementSpectrum` return the power and the tensor field. The parity fixture pins every
  polarisation, and its generator refuses an eigenvalue gap or a sense decided near its tolerance.
- Known answers: both sites of an element displaced 0.1 Å along the direct-lattice [1 0 1] of a
  monoclinic cell (14° from the Cartesian (1, 0, 1)) at q = (½, 0, 0) are found within 5°, with a
  share above 0.95; over 60 seeds the scatter is the sampling error σ/(δ√(mn)) rad per component.
  The existing planted modes are read along their axes.
- F3: the Šidák family is every element of the configuration with displacements, whichever element
  the call names. Asking about one element after seeing the data is selective inference, and
  naming it shrank the family to its own points: a weak doubling at 2.98σ over the 21 points of
  three elements read 3.30σ, "ordered", over its own 7. The payload's `family` names its
  `elements`. A call without an element is unchanged.
- `analyze_modulations` shows each ordered point's polarisation (a `polarisation` row: [u v w] to
  two decimals, share, the nearest <100>/<110>/<111> within 15°) and says it in the point's check.
  F2: a ratio or σ below its bar (3, 3σ) is never printed rounded onto it. Where the usual rounding
  reads at or above the bar, it is shown with the fewest significant digits that stay below (2.98,
  not 3), and nearer than six digits as `'<3'`; the no-order sentence adds "below the 3σ bar". The
  note says the correction is over every element's points. The tool gets its own 4,000-character
  result budget (`MODULATION_RESULT_BUDGET`). What the model sees changes in this tool's result
  only: the tool-result snapshot and a new sample of a call naming one element.
- Arm A+ reads a modulation claim's direction from the engine's polarisation, no longer from the
  PC1 of the element's most prolate site. On the dev set every arm A and D record, and every
  score, is unchanged; two A+ records move in their claim's direction and reading (re-recorded in
  `validation/baseline/arm-a-dev.json`). Approved by the manager on the user's delegation,
  2026-10-19 (`validation/findings/decisions/2026-10-19-mode-polarisation.json`).

### Copilot: one verification module answers "is this answer verified?" (architecture candidate #4)

- `src/llm/copilot/verification/` holds the evidence gate, the checks on the answer, numbers by
  reference and the computed Outcome behind one factory and two operations:
  `createVerification({ context, catalog, offered, phrasings, maxSendBacks, method })`,
  `record(step, show)` once per tool step (refs r1, r2, … in record order, none for a step with an
  error, whose facts are dropped; `show(ref)` builds the tool message and the ref resolves to
  exactly its data) and `judge(text, { canSendBack, rounds, stoppedBy, cut })` → a send-back, the
  turn's `done` event, or null for a blank reply with no answer standing. judge's order is fixed
  inside: fill the references, decide a send-back on the filled text, run the checks, add the
  unfilled-reference and missing-evidence findings, compute the Outcome. One verification serves
  one turn (record or judge after `done` throws); a phrasing for an unknown claim, or an asserted
  pattern that is not a global RegExp, is a TypeError.
- The internals are two declarative catalogues: CLAIMS (`claims.js`: one row per claim in
  send-back order — what, scope, phrasings, `backedBy` through `anyOf` / `firstOffered` /
  `allOf`, the prompt's phrase) and FACTS (`facts.js`: one row per facts key in check order — the
  producing tools, the shape, the phrasings it reads, the check; the `lowering` row has none yet).
  `evidence.js` and `answerAudit.js` are gone; `verdict.js` (with `isAnalysis`, the one
  definition of an analysis), `outcome.js`, `contextAuditFacts` and `shownContext` moved inside;
  `references.js` stays (without its unused `add`); the shell catalogue gives the finite clusters
  (`shells.js` `clusterElementsOf`, `CLUSTER_CN_MAX`). The turn builds one verification where it
  built four loop parameters, the method riding on it (both send-back wordings byte for byte);
  the agent loop records every tool step and judges every answer with it, and lost its reference
  store, step log, facts, gate state and `finish`.
- `EVIDENCE_RULES` stays written by hand; a binding test renders it from CLAIMS byte for byte
  (`evidenceSentence`) and pins what it does not say (`promptGaps`: coordination also takes
  `compute_bond_angles`, a lowering falls back to `detect_symmetry`, "no order" needs a per-site
  map and each finite cluster but neither `find_gr_peaks` nor `measure_symmetry_lowering`).
- Three behaviour changes, each pinned: the dock footer counts `tool_end.analysis` (views,
  unknown or skipped calls and a step stopped mid-run no longer count), as the computed Outcome
  does; the app's claim catalogue no longer reads the validation verdict's words — the harness
  passes them (`validation/verdict.js` `VERDICT_PHRASINGS`, the turn's `config.phrasings`), so
  its runs read as before; `DEFAULT_LIMITS` loses `maxEvidenceGates`, now the verification's
  `DEFAULT_MAX_SEND_BACKS` (2) — `config.limits.maxEvidenceGates` still sets it and the run policy
  still states it (`runPolicy.js` `TURN_LIMITS`, records unchanged). Also fixed: the dock fills an
  earlier segment's `{{ctx.…}}` from the context the model read.
- Unchanged: every send-back text, finding, computed verdict and `done` field; what the model
  sees (every model-visible snapshot); arms A and A+, the harness's expectations, the scripted
  replays and the answer corpus pins (its adapter moved onto the interface). No protocol bump.
  Tests at the interface (`verification.test.js`, `answerChecks.test.js`), the binding test
  (`verificationCatalogue.test.js`), the real tools' facts against FACTS on the Demo
  (`verificationFacts.test.js`) and the verdict's words (`validation/__tests__/verdictPhrasings.test.js`).
  Details: [the llm README](../web_app/frontend/src/llm/README.md), *Evidence, numbers by
  reference, the computed verdict*.

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

### Copilot validation: arm C without the method (plan item 1.8; protocol 2026-10-19.10)

- The turn takes `config.method: 'full' | 'off'` (`src/llm/copilot/method.js`, `TURN_METHODS`), and
  `planTurn` alone decides from it what the model reads. `'full'` (the dock, arm B) sends every
  request byte for byte as before. `'off'` (arm C, decided 2026-10-19) sends none of the method: no
  `read_skill` and no skills' index, no `EVIDENCE_RULES`, no method sentence of `TOOL_RULES` (the
  calibrated-σ reading, which analysis tells a distortion, effect size with z) or `ASK_RULES` ("Settle
  what the data can first: e.g. run find_gr_peaks …"), no method advice in a tool description or a
  check. The tools, their refusals, the Workbench's questions, the evidence gate and the answer checks
  run under both; `start.method` and the harness record say which one a turn followed.
- Tool specs gain an optional `method` field: the advice split out of the description, which keeps
  what the tool does and returns; put back with `'full'` (a string after a space, `{ after, text }`
  inside a sentence), left out with `'off'`. Split out: `list_sites`, `compute_bond_angles`,
  `analyze_clusters` (two pieces), `analyze_modulations`, `detect_symmetry`,
  `measure_symmetry_lowering`, `find_gr_peaks`, `compare_fit`.
- Checks built with `methodCheck` list the advice inside their text; arm C reads the fact and the
  verdict without it, and a refusal still refuses: the cut-window refusal, no angles, the pooled
  directions note, `analyze_clusters`' window failures and split-site note, `find_gr_peaks`' two
  split-shell notes; and three reading hints that name no tool (the manager, 2026-10-19):
  `find_gr_peaks`' "Ask whether the crystal has A–B distances there", `compute_bond_angles`' "Read the
  shape from histogram_10deg_sin_corrected …", `detect_symmetry`'s "Take inversion, chirality and
  polarity from this, not from memory". The text everyone else reads (the step cards, the record,
  arm A+) is unchanged.
- The shared system prompt's ground rules are written part by part (`src/llm/prompts/system.js`
  `GROUND_RULES`), each part a definition (what a context field is: quantity, units, computation,
  limits as a measurement), an inference (which field to use for a conclusion, or what a pattern is
  evidence of) or conduct (how to answer); a sentence that does both is split. Method `'off'` sends no
  inference part, the tool-free prompt included (the manager, 2026-10-19): of 26 sentences, 7
  definitions, 2 inference, 4 conduct, 13 split. `'full'` is byte for byte the prompt as before.
- The harness's `ARMS.C` sets `method: 'off'`; `ARMS.B` is unchanged. What the model sees moves for
  arm C only: `harness/first-request.arm-C.txt` (its system prompt is the `tools-ask.method-off`
  variant, 5,975 characters against 8,575), new `harness/tools.arm-C.json`, five
  `system-prompt/*.method-off.txt` variants and five `tool-results/*.method-off.json` samples; the
  dock's and arm B's snapshots and the budget pins pass unchanged. A test checks arm C's first request
  for every skill name, the skills' step headings, the method sentences, `EVIDENCE_RULES` and every
  piece of tool advice. Arms A, A+ and D are unaffected.
- The evidence gate checks the same rules under method `'off'`, but its send-back names no analysis
  and no tool: "… but this conclusion needs its evidence from an analysis run in this question. …"
  (`OFF_NUDGE`, decided by the manager on the user's delegation, 2026-10-19); arm B's wording is
  unchanged. The gate's notice and finding (the step card, the record) still name the analyses.
  Snapshots: `harness/gate-send-back.arm-B.txt` and `.arm-C.txt`.

### Validation: arm A+ reads the app's rules (architecture candidate 8; not a protocol change)

- Arm A+ (the engine sweep, `validation/armA.js`) applied three rules of its own where the app has
  one. It now calls the app's functions, held in `SWEEP_RULES`, so the control arm and the Copilot
  read a run the same way:
  - **The split site.** A+ used arm A's rule (PC1/PC2 ≥ 1.3 and κ₁ ≤ −0.5). It now uses
    `isSplitSite` (PC1 resolved, κ₁ ≤ −0.5 and at least 3 sampling errors below zero) with
    `list_sites`' static displacement δ = √(rms1² − rms⊥²) ± its error (`splitDisplacement`). Its
    readings say `list_sites` with κ₁, its z and δ. On synthetic prolate Gaussian sites, arm A's
    rule called 3.6 % (216 copies) to 14 % (100 copies) split; `isSplitSite` called none.
  - **The cluster elements.** A+ tried every same-element first shell peaking below 3.6 Å. It now
    tries the evidence gate's elements: those whose own first shell holds 0.5–4.5 neighbours. The
    rule is one function, `clusterElementsOf` (`llm/context/shells.js`, which now holds
    `CLUSTER_CN_MAX`), and `evidence.js` imports it.
  - **The lower group's tolerance.** A+ read the space group at the sampling noise floor. It now
    reads it at 1.25 × the floor: `lowerGroupTolerance` (`llm/context/runContext.js`), which is
    `measure_symmetry_lowering`'s default. The margin is `NOISE_FLOOR_MARGIN`, one constant read by
    the tool, the Workbench's tolerance question (`clarify.js`, until now its own `NOISE_MARGIN`)
    and A+.
- Arm A's rules are the frozen protocol's and are unchanged. On the dev set every arm A and D
  record is identical, and so is every score and claim. Every A+ record moves in its readings; on
  the ten paired cases:
  - The group is now read at 1.25 × the floor, with the same group in every case.
  - Eight cluster readings of extended networks (CN 6–12) are gone: the O of four perovskites,
    FeSi's Fe and Si, rock-salt Mg and O.
  - In dev-jahn-teller, two F sites are no longer split (2.8 sampling errors). The F claim stays.
  - dev-offcentre's δ is now 0.202 Å, where it was 0.201 (planted 0.2).
  - The judgment subset's five A+ records move in the same ways, with unchanged scores (their
    readings stay in the records and pins).
  - The moved records are re-recorded in `arm-a-dev.json`.
- What the model sees is unchanged: the model-visible snapshots pass as they were. Live runs stay
  comparable, and `PROTOCOL_VERSION` stays 2026-10-19.9.
- The skill lint gained a `rule` check (`RULES`). It flags a skill sentence that states one of the
  three rules otherwise. Its allow-list now holds five sentences, awaiting the user's skill review:
  - displacement-anisotropy SKILL.md:27 and significance-reading.md:17: a split site is a negative
    Mardia kurtosis;
  - coordination-shells SKILL.md:51: any polyhedron distortion goes to `analyze_clusters`;
  - local-symmetry SKILL.md:38: the floor itself is the line;
  - local-symmetry SKILL.md:59: the margin is left to the user.
- Known-answer tests: `validation/__tests__/aplusRules.test.js`. Decision record:
  `validation/findings/decisions/2026-10-19-aplus-app-rules.json`.

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

### Copilot: one shell reading per pair, per run (architecture candidate 6; protocol 2026-10-19.8)

- `src/llm/context/shells.js` reads the first coordination shell of every pair in the partials file
  once (`shellCatalogue(run)`, memoised per partials data, structure and run settings; frozen
  Shells), and every reader takes it from there: the run context's `pair_correlations`,
  `find_gr_peaks`, the Workbench's questions (`buildShellCatalog` is gone), the evidence gate, the
  answer checks (`cnMatches`), `compute_bond_angles`' window refusal (`countsCut`) and arm A+'s
  cluster windows. Before, five readers searched their own r ranges, rounded their own way and
  integrated the coordination number their own way.
- One rule per quantity: peaks below 8 Å and the closing minimum to max(8, r₁ + 3) Å, whatever range a
  caller passes; windows rounded to 0.01 Å and every CN integrated over exactly the window quoted
  with it (a count can move in its third digit: the Demo's Ta–Ta 2.76 → 2.75); separation and the
  closest-approach pile-up on one 2 % fraction (`find_gr_peaks`' FWHM edge clause is dropped; it
  moved nothing on the Demo or the dev cases); half an atom for every CN match.
- **The first peak alone ends at the minimum between the first two resolved peaks**, and exists only
  with two or more; r ± FWHM reached into the second peak. The Demo's Ta–Se first peak: 2.21–2.65 Å,
  3.25 Se per Ta (it read 2.33–2.71 Å, 3.86–3.9). Known answer: a 3 + 3 split with overlapping widths
  gives 3 ± 0.1 on three grids, where r ± FWHM gives 3.4–4.0. Approved by the manager on the user's
  delegation, 2026-10-19.
- `find_gr_peaks`: `first_peak_window_A` leaves; `first_shell.first_peak_only` (window, CN, trough)
  for every shell of two or more resolved peaks, its CN in the split-shell note; the shell no longer
  depends on `r_max` (a first peak above it gets a note and the shell, no longer "no peak"), as its
  description now says. The context adds `first_peak_to_A` / `first_peak_cn` for the same shells and
  lists peaks to 8 Å.
- The window-cut rule is unchanged: the Workbench asks on the CN margin, `compute_bond_angles` refuses
  on the coordination counts (60 % modal share, more than half unbonded; the Demo's 47 % text
  verbatim). Both readings live in the module; `validation/run.js shell-cuts` tabulates them on the
  dev windows (149: agree on 134, the counts alone refuse 11, the margin alone asks on 4). The
  CN-margin rule for the tool is a recommendation awaiting the user.
- The bond-angle engines (`triplets.py`, `workers/triplets.js`) return `coordinationSummary`
  `{mean, modal, modalShare, unbondedShare}`, the lowest count winning a tie; the Bond Geometry page
  and `compute_bond_angles` read it instead of recomputing it. Both golden fixtures pin it.
- The skill lint gains a value check (Demo values a skill quotes, against the tool): the four places
  the coordination-shells skill quotes the old first peak are allow-listed with both values until the
  user's skill review, a precondition for the freeze. What the model sees moves (`find_gr_peaks`, the
  context, the tool's description); the `compute_bond_angles` snapshots do not; the dock's worst
  turn grows by 90 tokens, arm B's by 16.
  Arms A and A+ on dev: scores and claims unchanged, three A+ cluster windows 0.01 Å lower at their
  start (re-recorded).

### Tests: an App-level smoke test of the Copilot dock's state lifecycle; the Demo no longer follows a folder's Live Data (plan 1.4)

- `src/__tests__/appDockLifecycle.test.jsx` drives the real App, Dashboard, Detected SG card, run
  module and dock on the Demo run's files, with no model (a recording fake at
  `streamCopilotChat`, where the dock's model port sends each request) and no server (axios refuses
  every call; fetch serves only the Demo's files and a stubbed model list). The Demo's files are
  picked as a folder in the static build with Live Data on; RMCProfile's save is a new `.rmc6f` in
  that folder, seen by the Live Data poll.
- It checks, end to end: after a save the dock draws the configuration divider; the next request
  carries no earlier question or answer, and its system prompt, tools and context are a fresh
  conversation's (the context is the run as App hands it now, the first context with only the new
  move count); memory resumes from the divider; a Detected SG pick (a ladder brick) is kept across
  another save of the same run, on the card and in the Copilot's context; another run (the Demo
  button) brings the card and the context back to `DEFAULT_VIEW.toleranceA` and a new conversation.
- **Fix: opening the Demo stops static-mode Live Data watching the folder picked before.** Before,
  a later save into that folder replaced the Demo's data (name, plots, structure) under the Demo's
  `runKey`, so the Copilot kept the conversation and the Detected SG pick. `handleToggleDemo` now
  clears the folder handle and its signature, and each poll tick reads the folder open at that
  moment and drops a read whose folder was replaced meanwhile. Live Data switched on with no
  folder open now watches the next folder picked (before, it did not start until switched off and
  on). Pinned by two more cases in the same test file.
- Not covered in jsdom (left for the browser check): the real File System Access handle and Web
  Worker (stand-ins here), the provider's streaming, the dock's sheet layout, and Flask-mode Live
  Data with the real dock.

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

### State lifecycle: a new configuration ends the Copilot's memory; the Detected SG pick resets with the run (plan 1.4)

- **A Live Data save that brings a new `.rmc6f` is a new configuration of the same run.** App
  passes the AI Copilot a `configKey` beside its `runKey`: `configEpoch` in Flask mode, and for a
  browser-loaded run (static-mode Live Data included) the signature of its `.rmc6f` files
  (`configurationSignature` in `browserData.js`, which the Flask check now uses too). A save of the
  plots or the `.log` alone changes neither.
- **The conversation stays, but the model no longer reads the old answers as current.** Before, a
  save kept up to four earlier exchanges about the previous configuration in the model's memory
  while the context described the new one. Now each question remembers the configuration it was
  asked under; the dock draws a one-line divider where it changed ("New configuration saved —
  earlier answers were about the previous one", announced to screen readers), after an answer that
  was still running at the save, and a turn's history holds only the questions asked under the
  current configuration since the last change (`turnHistory(turns, configKey)`,
  `configurationBreaks`). The choice keeps the user's conversation readable without misleading the
  model; starting a new conversation on every save would have thrown away what the user can still
  read.
- **What the model is sent changes in one case only:** the first question after a configuration
  change is sent with no memory, exactly as a fresh conversation's first request, and later
  questions remember only the exchanges since the change. The first request of a fresh
  conversation is unchanged (the model-visible snapshots pass as they were). The validation asks
  one question per fresh run, with no Live Data save, so its requests are unchanged: no protocol
  bump.
- **The Detected SG pick resets with the run.** The shared tolerance (`useRunTolerance` in
  `symTolContext.js`) returns to the card's default (`DEFAULT_VIEW.toleranceA`, 0.2 Å) whenever
  another run opens (`runKey`), even one that was open before, so one run's pick never sets the
  next run's card or Copilot context. It is kept across page switches and across a Live Data save or
  a Load of the same run. The analysis pages' Live Data behaviour is unchanged.
- Pinned by `src/__tests__/appStateLifecycle.test.jsx` (Flask and static-mode Live Data: run vs
  configuration, the pick kept or reset), the dock's divider and memory in
  `src/llm/__tests__/copilotChat.test.jsx`, `turnHistory` / `configurationBreaks` in
  `src/llm/__tests__/turn.test.js`, and `configurationSignature` in `browserData.test.js`.

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
  cluster and modulation readings change (the ten dev A+ records, re-recorded in
  `validation/baseline/arm-a-dev.json`). Only a held-out case would exercise the rule changes.
  Approved by the manager on the user's delegation, 2026-10-19.

### Copilot: the provider quirks behind the model port (architecture candidate 3)

- `providerModel` (`src/llm/provider/model.js`) owns what differs between providers; the turn's
  and the loop's interfaces are unchanged. Its `prepare` is now async: it detects the server
  (`detectServer`, which probes only localhost or port 11434), puts the tool choice in that
  server's form, and asks the context guard for its decision (`planCopilotChat`, split out of
  `streamCopilotChat` in `contextGuard.js`). It returns the request as it goes out, with
  `plan`: the server kind; Ollama's `numCtx`, `numPredict`, estimate, loaded and maximum context;
  LM Studio's loaded context; and the tool choice's form (`TOOL_CHOICE_FORMS`: `body`, `omitted`,
  `tools_withheld`). `stream` sends under that plan. The guard's numbers and rules are unchanged
  (2.5 and 1.4 characters a token, the 8,192-token reserve, LM Studio's 4,096 tokens of room).
- **Behaviour change for Ollama runs:** the final request (the analysis budget used, `tool_choice:
  'none'`) now carries no tools, because Ollama's native `/api/chat` has no `tool_choice`. The
  model can no longer answer it with calls only, so it no longer ends on the "budget used" fallback
  for that reason. The earlier calls and results still travel as messages. The estimate counts
  what is sent, so the final request's estimate leaves out the tool schemas. OpenAI, Gemini,
  Anthropic and LM Studio get the same body as before: `"tool_choice"` beside the tools.
- **Typed failures.** A request the guard refuses throws `ContextLimitError` from `prepare`, before
  anything is sent. The turn then emits no `request` event for it, as for a spend cap, and a record
  has no entry for it. A model that cannot take tools fails `stream` with a
  `ToolsUnsupportedError`, read from the provider's response: Ollama's 400 `{ error: "<model> does
  not support tools" }`, OpenAI's 400 `{ error: { param: 'tools', code: 'unsupported_parameter' } }`,
  then the old text pattern as the port's one fallback (LM Studio, Gemini, Anthropic). The agent
  loop retries on the type, not on the text: `isToolsUnsupportedError` and its regex left
  `agentLoop.js`. `client.js` and `ollama.js` keep the parsed error body on the Error (`httpError`).
- The loop keeps its final-round guard: calls in the final round are never run, and a reply with no
  answer ends on the fixed "budget used" answer. LM Studio, Gemini and Anthropic take the field but
  are not verified to honour it. The port's server probes now run in `prepare`, which the loop
  puts under the stall clock (`beforeStall`, `robustness.js`).
- Validation (protocol 2026-10-19.7): the run policy's Ollama final-round form is now "no tools
  sent", and the scripted check pins it. Each request's `num_ctx` is read from `sent.plan`
  (`observeContext` is removed). Ollama runs under earlier protocols are not comparable on their
  final round.
- Tests: `modelPort.test.js` (the port over fake servers: the plan and form per server, the refusal
  in `prepare`, the typed failure and every case it must not catch, the loop over the port);
  `planCopilotChat` in `contextGuard.test.js`; the refused request with no `request` event in
  `turn.test.js` and `runPolicy.test.js`. The model-visible snapshots, `arms.test.js` and the
  budget test are unchanged. Details: [the llm README](../web_app/frontend/src/llm/README.md), *The
  model port*, and [copilot-validation.md](algorithms/copilot-validation.md), Protocol 2026-10-19.7.

### Copilot validation: a fit-data dev case (plan 4.6)

- **dev-fit-unconverged**, a dev-only case (`devOnly` in `validation/cases.js`; no test twin, no
  reported number), is the first case with fit files, so `assess_convergence`, `compare_fit` and
  the rmc-fit-assessment skill are exercised. Rock-salt MgO, 6³ cells, isotropic thermal motion
  only (no order, Fm-3m), seeded. Besides the `.rmc6f` (now with RMCProfile's run counters: 9.4
  accepted moves per atom) and the partials, its run folder holds RMCProfile's X-ray file family
  in the Demo's formats: `rmc_run_FQ1.csv`, its Fourier transform `rmc_run_FT_XFQ1.csv` and the
  χ² log `rmc_run-00.log`. `runFolder.js` passes them as run-folder entries (`caseFiles`), so the
  run module detects their kinds as for a user's run.
- The planted fit truth (`truth.fit`): the X_ray_(R)1 χ² is still falling after a fast drop (the
  plot's end looks flat; the watchdog and `assess_convergence` say improving, ln χ² −0.034 over
  the last tenth), and the measured G(r) holds a peak at 2.55 Å the configuration does not
  reproduce (`compare_fit`: 97 % of the G(r) misfit in 2.17–3.33 Å; in F(Q), the same
  measurement, a misfit at low Q).
- Scored on its structure like every case; the scoring and the verdict schema are unchanged. The
  fit truth is checked by the context, the fit tools and the answer check's convergence rule
  (`validation/__tests__/fitCase.test.js`, scripted model). Arms A and A+ score 1 on it (appended
  to `validation/baseline/arm-a-dev.json`); the other 27 dev records and the test set are
  unchanged.
- Pins updated: `arms.test.js`, `split.test.js` (a dev-only kind), `caseBlind.test.js` (the dev
  cases now offer every tool its sweep calls), `toolDataShapes.test.js` (the fit tools on the
  case), and `contextBudget.test.js`: arm B's largest dev case is now this one, 35,906 tokens
  needed, 3,138 over 32,768.

### Copilot validation: the run policy, and LM Studio models at 48k (plan 1.2)

- `validation/runPolicy.js` states what a live run runs under besides the cases, the prompt, the
  tools and the scoring: the server kind and version, the temperature and the warmer retry's
  floor, the loop's limits, the minimum window per server kind, the window read before and after
  the run, the `num_ctx` each Ollama request was sent, the reply cap, the context guard's
  decisions and the final round's `tool_choice` form. Each record keeps it (`record.policy`, and
  `num_ctx` on each request), each transcript prints it under its header, and the report adds
  *What a comparison across models also compares*: the columns' policies side by side, whether
  the limits, temperature and guard are the same in every run, and what no policy sets (sampling
  beyond the temperature, the reasoning mode, the tokenizer, the quantization).
- The decision of 2026-10-19 (the manager, on the user's delegation): the guard, its estimate, the
  8,192-token reply reserve, the budgets and the prompt are unchanged, since refusing rather than
  cutting is wanted. The live runs load each LM Studio model with at least 49,152 tokens (48k),
  which holds the worst two-skill turn of the dock and of arm B with 10,590 and 14,033 tokens to
  spare (`contextBudget.test.js`). Before every run, `run.js live` reads the server
  (`preflight`): another kind, none, a model not loaded, or one loaded below 48k stops the batch
  with what was found, what is needed and how to reload (`lms load <model> --context-length
  49152`). Nothing is sent and no record is written, so the batch resumes. Ollama still sizes
  `num_ctx` itself and is recorded only. `models.js` names each model's server kind.
- A scripted check sends a final round to each provider preset through the production port and
  pins what goes out today: `"tool_choice": "none"` with the tools for OpenAI, Gemini, Anthropic
  and LM Studio, and no `tool_choice` on Ollama's native API, where the tools are still offered.
  No provider behaviour changed; enforcing the final round on Ollama is the model port's job.
- Protocol 2026-10-19.5: what the model is sent is unchanged, but a run's conditions are now part
  of the protocol and earlier records do not state theirs. Details:
  [copilot-validation.md](algorithms/copilot-validation.md), *The run policy*.

### Copilot: one module runs one question; the dock and the harness are its adapters

- `src/llm/copilot/turn.js` is the turn module: `runTurn(spec)` takes the question, the run props,
  a model port, the caller's tools, an ask-user port, the skills, the history, a signal, a clock
  and `{ toolsEnabled, method, limits }`, and builds everything else — the context and the context
  as the model reads it, `ask_user` and `read_skill` after the caller's tools, the prompt variant,
  the history window, the Workbench's questions, the evidence gate, the answer checks' context
  facts and the references. Its events are the loop's, each stamped with `at`, between a `start`
  (what is offered, the prompt, the method, whether tools are active, the shown context) and
  exactly one terminal event: `done`, `stopped` or `failed` (with the message and the hint).
  `describeTurn(spec)` is the same assembly, not run. Invalid input (a tool named `ask_user` or
  `read_skill`, no model) throws at the first step; after `start` it never throws.
- `src/llm/provider/model.js` is the model port (`prepare(intent) → sent`, `stream(sent)`), and
  `providerModel` its production adapter over `streamCopilotChat`: a request asked again after a
  looping reply runs at a temperature of 0.6 or more, and the request as sent never holds the API
  key. The agent loop calls the port and emits each request as sent, with why it was made.
- The dock (`useCopilot.js`) keeps its React state, Stop, the question card, the card focus, the
  models that refused tools and the run reset; CopilotChat no longer builds the clarifier, parses
  the skills or adds `ask_user` and `read_skill`. What the dock sends is unchanged (the
  model-visible snapshots pass as they were). `CHAT_HISTORY_TURNS` is renamed `HISTORY_MESSAGES`.
- The validation harness runs each question through the same module (protocol 2026-10-19.3), with
  four changes, each pinned: the warmer retry at 0.6 or more (it ran at the model's 0.2), the
  no-tools prompt for a model that refuses tools, number provenance against the context as the
  model read it, and every request kept in the record (`record.requests`). `run.js` builds the
  port with `providerModel`, and `spend.js` `meteredModel` meters it (replacing `meteredCall`; a
  request the cap refuses is never sent). Arms A and A+ are unchanged.
- Tests: `turn.test.js` at the interface (a scripted port, a fake ask-user, a fake clock);
  `modelVisible/turnRequests.test.js` (the turn's first request for a dock-shaped and a
  harness-shaped spec equals the snapshot files); the prompt-assembly, memory, warmer-retry and
  "when to ask" tests moved onto the turn; the scripted replays run through `runTurn` with their
  assertions unchanged. Details: [the llm README](../web_app/frontend/src/llm/README.md), *The
  turn*, and [copilot-validation.md](algorithms/copilot-validation.md), Protocol 2026-10-19.3.

### Copilot: a corpus of real answers for the gate's and the checks' claim readings (plan 3.5)

- `src/__tests__/answerCorpus.test.js` runs real model answers through the evidence gate and the
  checks on the answer as they are today. The answers sit one file each in
  `src/__tests__/answerCorpus/entries/`, each with its text, source (file, lines, commit), model,
  question and run (the Demo or a dev case, never a test case). Where the record shows them, an
  entry also lists the steps that ran, and those steps must reproduce the transcript's numbers
  with the real engines.
- One adapter (`answerCorpus/readAnswer.js` `readAnswer`, public functions only) gives the reading
  each entry pins: the claims the gate reads, its send-back given the steps (or "steps
  unknown"), and the checks that fire. Each reading carries a person's judgement, so false alarms
  are counted.
- Nine answers today: gemma4:26b's three correct ones of 2026-10-11 (verbatim), Muse Glimmer's
  GaSe₄ answer from the README GIF (transcribed), and five wrong ones of the 2026-10-10 spot check
  and review, as the docs and skills quote them. Four readings fire, all right: chirality,
  inversion lost, convergence, and the gate's cluster send-back. There is no false alarm and one
  miss ("no cluster distortion" from a pooled map needs no analysis). No regex, check, prompt or
  snapshot changed.
- [copilot-validation.md](algorithms/copilot-validation.md#the-answer-corpus-plan-item-35)
  says how to add an answer from the live runs' records. A regex change must keep the pins or
  update them with a judgement.

### Copilot: `list_sites` shows the split (plan 4.1)

- Each `list_sites` row now carries, where it applies, the excess kurtosis along each resolved
  principal axis (`kurtosis_axes`, null where an axis is not resolved), PC1 as a direct-lattice
  direction (`pc1_uvw`, when PC1 is resolved), and for a split site its static displacement with its
  sampling error (`delta_A`, `delta_se_A`), so an answer quotes the engine's number by reference
  instead of working out √(rms1² − rms⊥²) by hand. A site is split by the rule the cluster summary
  already used (PC1 resolved, PC1 excess kurtosis ≤ −0.5 and three sampling errors below zero), now
  one function, `isSplitSite`, for both tools.
- δ² = λ₁ − (λ₂ + λ₃)/2: rms⊥ is the RMS of the two other axes, unbiased where PC2 and PC3 are
  equal (rms2 alone, the larger of two equal sample variances, reads δ low). Its error comes from
  the eigenvalues' sampling errors, (κ + 2)·m₂²/n each, by the delta method
  ([ai-assistant.md](algorithms/ai-assistant.md), Step 20). Known answer over 1,000 seeds for each of
  five planted splits (equal and unequal occupation, 216–512 copies, anisotropic thermal motion;
  cubic, hexagonal and monoclinic cells): (δ̂ − δ)/SE has a mean of −0.01 to 0.07, an sd of
  0.99–1.04 and 94–96 % within ±1.96, and PC1 lies within 7.6° of the planted [u v w]. No Gaussian
  ellipsoid is called split. On the dev cases the planted 0.2 and 0.1 Å splits read
  0.202 ± 0.0028 and 0.0989 ± 0.0029 Å.
- Split sites come first in the table, so trimming to the result budget never drops one before
  another site, and a note check names them; the tool's description lists the new fields. What the
  model sees changes: the `list_sites` result and description snapshots. The description adds 38
  tokens to every request, so the 32k budget test's pins move: the dock's worst turn needs 38,562
  tokens (5,794 over 32,768) and arm B's 35,119 (2,351 over). Arms A and A+ do not read
  `list_sites`, and their dev records are unchanged. Approved by the manager on the user's
  delegation, 2026-10-19.

### Copilot validation: one file per decision record

- Decisions on the validation itself are kept one JSON file per record,
  `validation/findings/decisions/<id>.json` (the id: the date and a slug), read in id order by
  `parseFindingFiles` and validated by the findings schema, so parallel changes never edit the same
  file. `decisions.jsonl` is gone: its one record (the recomposed held-out breathing case) moved
  there unchanged, beside plan 4.1's, `2026-10-19-list-sites-split.json`. `run.js findings <folder>`
  validates such a folder.

### Copilot: the budget test for a 32k window (plan 5.5)

- `src/__tests__/contextBudget.test.js` runs the worst realistic turn through the real agent loop
  with a scripted model and counts it with the context guard's own estimator and constants: the
  dock's largest system prompt and full tool set, the Demo's context block, the two largest
  skills read whole (`read_skill`'s budget and keep rule), a round of analyses filling the loop's
  12,000-character turn budget, and the 8,192-token reply reserve. The same for the harness's arm
  B on the dev case with the largest first request.
- It does not fit today: the guard counts tool messages at 1.4 characters a token, so the dock's
  worst turn needs 38,525 tokens (5,757 over 32,768) and arm B's 35,081 (2,313 over). The test pins
  every term and the overflow, and keeps "fits in 32,768" as an expected failure (`it.fails`) that
  turns red once it fits. Plan 5.5's other half, moving long skill material to on-demand
  references, waits for the user's skill review; no prompt, tool, skill or snapshot changed.

### Tests: a lint for the Copilot's skills (AGENT_PLAN 5.1)

- `src/__tests__/skillLint.test.js` reads every `.claude/skills` file the app bundles (it never
  edits them) and checks three things. **Tools:** every snake_case name in a skill is a tool the
  builders offer on the Demo (all 16), or else a parameter or a result or context field the model
  sees. **Gate numbers:** each threshold a skill states equals the constant of the module that
  decides it (`skillLint/tables.js` `BINDINGS`: 16 bindings to `SIGNIFICANCE_SIGMA`,
  `MODAL_COORDINATION_SHARE`, `CLUSTER_Z`, `MAP_TEST_MIN_PAIRS`, `AXIS_RESOLUTION_SIGMAS`, …, some
  with the statistic they gate and a probe of the code), and every number after a comparator is
  bound or listed as no gate. **Firm rules:** each "*Firm rule in the app:*" is mapped to the
  answer check, refusal or context that enforces it, and that code is run.
- Today's three violations, all in displacement-anisotropy, sit in an allow-list with why: the
  low-index cone has no number (the app uses 15°), "preferred direction" does not name the peak
  test, and the split site is a negative kurtosis where the tools use PC1 ≤ −0.5 (plan 5.3). A new
  violation fails, and so does an entry a fix has cleared.
- `workers/pcaKde.js` exports `AXIS_RESOLUTION_SIGMAS` (no change in behaviour).

### Copilot: one module for the run as the Copilot sees it; the harness records its view

- `src/copilotRun.js` turns an open run and the view on screen into everything the AI Copilot is
  handed: CopilotChat's eight run props, the engine tools, the elements, the engine and the
  normalised view (`DEFAULT_VIEW`: the Detected SG card's 0.2 Å, nothing hidden, no PCA table,
  Live Data off). Its browser adapter `useCopilotRun.js` takes over the Dashboard's symmetry and
  `.dat` work and copilotTools.js's `useCopilotTools`; its in-process adapter `openRunFolder.js`
  (Node and tests only; an eslint rule keeps it out of `src/`) is how the validation harness now
  opens each case. The Dashboard's plot list, χ² history and plot parse are shared from
  browserData. App behaviour and what the model sees are unchanged (the model-visible snapshots
  pass as they were; arms A and A+ are unchanged).
- Plan item 1.3: each live record keeps the view it ran with (`record.view`) and the context
  sections the model had and lacked, with `earlier_turns: 0`; `run.js live --pca-table` (dev
  only, off by default) adds the PCA table for a check of its effect
  ([copilot-validation.md](algorithms/copilot-validation.md#records-transcripts-findings)).

### Copilot validation: the held-out breathing case recomposed

- The test breathing case's cluster and anion were the skills' reference compound, whose
  published shells `coordination-shells/lacunar-spinel-shells.md` gives and `read_skill` serves to
  arm B: a leak of the case's kind to arm B only. The case is recomposed, not the skill (the
  user's text): a composition a fresh author chose without dev results, which no skill (main or
  PR #71's drafts) or model-visible snapshot names. Its seed follows a rule fixed before looking:
  the first seed of a fixed sequence whose drawn breathing sd over the case's clusters lies in a
  fixed window about the planted sd. Shape, breathing, thermal width and supercell are unchanged.
  (Values, seed and rule sealed with the test cases since the holdout gate; see `TEST_SET_SHA256`.)
  Approved by the manager on the user's delegation, 2026-10-19; the user may reverse it.
- `TEST_SET_SHA256` and `PROTOCOL_VERSION` 2026-10-19.2. No arm B or C run had touched the test
  set. Arms A and A+ re-recorded on it with `--final` (4ab9eab): only the breathing case's records
  change, and no score moves (A 0.81, A+ 0.99); the dev records are unchanged.
- The holdout scan's open-leak list is empty (kept for future skill-text hits), and its four
  allow-list entries for the old cluster are gone. `validation/findings/decisions.jsonl` records
  the decision in the findings' form, validated by `split.test.js`. Details:
  [copilot-validation.md](algorithms/copilot-validation.md), Protocol 2026-10-19.2.

### Copilot: the shape of every tool's data, declared and checked (plan 3.3)

- With numbers by reference the model quotes `{{r3.sites.0.shift_A}}`, so the paths of a tool's
  `data` are an interface. Each live tool's spec now declares them: `dataShape`, in a small
  notation that reads like the JSON (`'number|null'`, `'key?'` optional, `'cn_*_per_*'` pattern
  keys, `[shape]` lists, `[a, b]` tuples, `oneOf(…)`), with one pure checker
  (`src/llm/copilot/dataShape.js`) that names every mismatch by its reference path. The four run
  tools (with every `read_run_data` section), the seven engine tools, `ask_user` and `read_skill`
  declare theirs; the engine tools' shapes sit beside the summarizers that build the data.
- `dataShape` is internal: `toRequestTools` sends name, description and parameters only, and the
  model-visible snapshots are unchanged.
- `src/__tests__/toolDataShapes.test.js` runs every tool through `executeToolCall` on the Demo run
  and on every dev case (≈ 3 s), with arguments read from each run's own data (every section and
  pair, whole, cut, confirmed, two-window and empty bond windows, each same-element shell's
  clusters, three tolerances, each `ask_user` reply form, every skill and reference), and checks
  each result both ways: no undeclared key, no required key missing, every value of its type, and
  every declared key returned by some call except the forms the runs cannot reach (listed with
  their reason). A PR that changes a tool's output changes its `dataShape` in the same diff. No
  tool's behaviour or output changed.

### Copilot validation: the holdout scan widened

- `validation/__tests__/split.test.js` (with `holdoutScan.js`) now looks for each held-out test
  case's composition, cluster, structure type, planted lengths, directions and wavevector, read
  from `cases.js` at run time, not only its id. It scans `src` (js, jsx, md and the model-visible
  snapshots), every skill file, `rmc_toolkits`, `docs/` and the top-level docs; the validation and
  copilot-validation.md are not scanned. A bare element symbol is never a hit, and a failure prints
  the matched text masked. Only text no model can read may be allow-listed, and the test checks
  that rule; five test- or doc-only hits found on 2026-10-19 are. The skills' reference compound
  reaches arm B through `read_skill`, so it is an open leak pending the user: printed on every run
  and pinned line by line. Both lists are in
  [copilot-validation.md](algorithms/copilot-validation.md#the-holdout-scan) (plan item 1.5).

### Copilot validation: the harness no longer tells the model the case kind

- The harness named each run, and its `.rmc6f` and partials files, after the case id, and the
  run name reaches the model in its context (`"run": "dev-offcentre"`): every arm B and C run was
  told its case kind. Every case is now the run `rmc_run` (`validation/runFolder.js` `RUN_NAME`,
  `RUN_FILES`); the record keeps the case id. `run.js write` saves each case's run folder as
  `<out>/<case>/rmc_run/`, with `truth.json` beside it rather than in it.
- `PROTOCOL_VERSION` 2026-10-19.1: live runs under an earlier protocol, if any, are not comparable.
  Arms A and A+ are unchanged. The model-visible snapshot changes in one line, the context's `run`.
- `src/__tests__/modelVisible/caseBlind.test.js`: no case id or kind in the harness's first
  request on any case (dev and test, arms B and C), and no case id in any tool result on the dev
  cases ([copilot-validation.md](algorithms/copilot-validation.md#what-the-model-sees-snapshot-tests)).

### Tests: generated goldens at the orientation and site seams

- The displacement-direction engine and its JS port are held to one generated golden instead of
  values copied by hand into both suites: `tests/generate_orientation_fixture.py` writes
  `orientation_fixture.json` (21 full payloads with their key sets — the golden cloud, the 8-well
  model, isotropic, one-sided, prolate in the PCA frame, sparse, singly occupied pairs and demo-run
  sites, at ν 2–10 — plus the tie-family cells, neighbour checksums, polygon starts, the
  `recommendedFrequency` table and Pearson's null moments), checked by `orientationParity.test.js`.
  The hand-copied tables (`GOLDEN_PEAK`/`_MAP`/`_ASYMMETRY`/`_ANISOTROPY`, `WELL_MODEL_PEAKS`,
  `SPECIAL_DIRECTION_CELLS`, `AXIS_CELLS_NU3`, `NEIGHBOR_CHECKSUMS`, `POLYGON_STARTS_NU38`,
  `PEARSON_MOMENTS`, the `recommendedFrequency` pins) are gone; the generated values equal them bit
  for bit.
- The PCA site table (`site_ellipsoids` ↔ `siteEllipsoids`) gets its first cross-engine golden:
  `tests/generate_pca_sites_fixture.py` → `pca_sites_fixture.json` (all 52 demo-run sites and
  synthetic isotropic, split, zero-spread and mixed sites), checked by `pcaSitesParity.test.js`.
- `tests/test_orientation_fixture.py` and `tests/test_pca_sites_fixture.py` re-run the generators
  on the committed inputs, so a fixture cannot go stale. Both sides compare by one rule
  (`goldenCompare.js` ↔ `tests/golden_compare.py`): integers, booleans and nulls exactly, other
  numbers to 1e-9 relative + 1e-12, p-values to 1e-9 relative, eigenvector sets eigenspace by
  eigenspace (a degenerate one, where LAPACK and the JS Jacobi solver pick different bases, by its
  projector). No engine changed.

### Copilot: snapshots of what the model sees

- `web_app/frontend/src/__tests__/modelVisible/` pins, as readable files, every system-prompt
  variant, the first request and tool list of the dock (Demo run) and of the validation harness
  (`dev-offcentre`, arms B and C), and one tool message per tool on the Demo. Tests only; a diff in
  them is a change to what the model sees, to be called out in the PR
  ([copilot-validation.md](algorithms/copilot-validation.md#what-the-model-sees-snapshot-tests)).

### Copilot: no conclusion without its evidence; numbers by reference; the verdict computed (THY-18)

- **Evidence gate:** an answer that concludes order between cells, a split site as disorder, a
  preferred direction, a cluster distortion, a coordination number, a lowering or "no order" must
  have run, in that question, the analysis that shows it; otherwise the agent loop withdraws it and
  sends it back once per rule (at most twice a turn) with what to run. It is a check after the
  answer, since not every model takes a forced tool choice; evidence still missing is a finding.
- **Numbers by reference:** every tool result carries a `ref`; the model writes
  `{{r2.peak.significance_sigma}}` and the app fills in the engine's value where the answer is
  shown, recorded and scored. A reference that names nothing is a finding.
- **Computed verdict:** the Outcome chip comes from the checks (not verified, not achieved,
  partial, achieved); the model's Outcome line only explains, and its own grade is kept as
  `claimed`. `analyze_clusters` gets a 4,500-character result budget so its numbers reach the model
  whole.
- The validation runs with all three (protocol 2026-10-17.2) and reports the gate's send-backs, the
  references and the computed verdict against the score beside the case score.

### Copilot validation batch 1: the review's tool changes

Every tool change the 2026-10-14 tool and skill review proposed ([copilot-review.md](algorithms/copilot-review.md)),
approved by the user that day, each with a regression test and a finding in
`web_app/frontend/validation/findings/batch-1.jsonl`:

- **L1, split shells:** a peak within 0.5 Å of the first joins a shell cut at a deep trough when the
  summed CN is 4, 6, 8 or 12 (a Jahn–Teller CuF₆ is CN 6, never four-coordinate); the first peak
  alone is offered and asked about.
- **L2, the cluster null:** `analyze_clusters` reads the shape at each cluster's own breathing (a new
  shape null in `clusters.py` and its port, with the radial noise the per-cluster fit takes out put
  back) and the breathing against independent atoms; split member sites are named.
- **L3, directions:** a preferred direction needs the peak test; the map test is anisotropy.
- **L4, two analyses:** `measure_symmetry_lowering` gives the coherent shift of a lowering (Å,
  direction, against the noise floor); `analyze_modulations` also tests general wavevectors.
- **L5, answer checks:** cluster claims (distorted, breathing, rigid) and coordination numbers.
- **Skills** (a PR of their own, for the user's review): the four drafts take the review's
  proposals, and a new draft, `local-order`, says when "no order" may be concluded. All five await
  the user's interview.

The A+ control follows the new engines (dev-symmetry-lowering 0.83 → 1).

### Charts on touch: a finger scrolls the page; a sideways swipe zooms x, a tap reads the value

The chart SVG carried `touch-action: none`, so on phones and tablets every finger drag on a chart
box-zoomed it, and a phone Dashboard, which is mostly charts, could hardly be scrolled. It is now
`touch-action: pan-y pinch-zoom`: a vertical drag scrolls the page and a pinch zooms the page, as over
any other content. A drag that starts sideways selects an x range and zooms x only, with y re-fitting
to it (a finger drifts far more than the 8-unit box threshold, about 4 CSS px on a phone card, so a
box would squash y on nearly every swipe). A press that travels at most 16 CSS px is a tap: it shows
the value readout, which touch never had, until the next press; a double tap resets. The readout's
crosshair is marked `data-export="omit"` and left out of saved figures. The mouse keeps box zoom,
hover, double-click and ctrl/⌘ + wheel (`interactivePlotTouch.test.jsx`; checked with CDP touch
input in headless Chromium at 375 × 812 on the Dashboard and Bond Geometry).

### Charts: the wheel scrolls the page; ctrl/⌘ + wheel or a pinch zooms

React registers `onWheel` as a passive listener, so a chart's `preventDefault()` was ignored
(Chromium logged "Unable to preventDefault inside passive event listener invocation" on every wheel
tick) and a wheel over any chart both scrolled the page and zoomed the chart. A plain wheel now only
scrolls the page; ctrl/⌘ + wheel, and a trackpad pinch (which browsers send as ctrl + wheel), zoom x
through a native non-passive listener without zooming the page. The zoom is anchored at the pointer,
a notch is 1.22× either way (exact inverses), and a pinch's small deltas zoom in proportion and
compose between renders. Drag-to-zoom and double-click are unchanged
(`interactivePlotWheel.test.jsx`; checked with trusted wheel events in headless Chromium on the
Dashboard and Bond Geometry).

### README: the visual front page back, with a Copilot GIF

The one-screen README (2026-10-11) kept the facts but lost what made the page inviting. It gets
back the three-step *Open the app* block, the 🔒/⚡/🤖/📖 callouts and the screenshot grid, and
leads with a new GIF of the AI Copilot answering a question on the Demo run with a local model
(`assets/rmc-toolkits-copilot.gif`, which replaces the chat-only `rmc-toolkits-assistant.gif` in
QuickStart too). The feature essays, the tested-models table and the signature equations stay in
the Copilot guide and ALGORITHMS.md.

### AI Copilot: robust on local models

For the local models the Copilot runs on (qwen3.8-27b on LM Studio, a dense reasoning model that
took 1–5 minutes a question in the 2026-10-10 spot check; gemma4:26b on Ollama), a turn no longer
hangs, loops or ends in nothing without saying why (`llm/copilot/robustness.js`):

- **Deadlines.** A model request that sends nothing for 3 minutes (`stallMs`; reasoning counts as
  progress), or runs past 10 (`requestMaxMs`), is cancelled, and the error says what to do.
- **Loops.** A reply or its reasoning that writes the same passage again and again is stopped as
  soon as it is caught and asked again once, at temperature 0.6 or more; two rounds in a row whose
  calls only repeat earlier ones or fail end the analyses (`stoppedBy: 'no_progress'`), so
  repetition never spends the round budget.
- **Empty replies.** A reply with neither an answer nor a call (a reasoning-only reply) gets one
  reminder before the turn ends, and the footer says why it ended empty. A call written as JSON
  text is run as a call. Notes to the model join the last user message, as Gemma's chat template
  refuses two user messages in a row.
- **Compaction keeps the numbers.** An earlier result shortened to fit the turn keeps its top-level
  numbers next to its checks, and is never shortened twice; an earlier answer re-sent as memory
  is cut to its sentences with numbers and its Outcome line above 2,000 characters.
- **Errors explained.** Out of memory, a context overflow, no model loaded and a lost connection
  each get a line on what to do.

### Copilot validation review (2026-10-14)

A review of the Copilot's tools against synthetic configurations with planted local structure
(`web_app/frontend/validation/`, [docs/algorithms/copilot-validation.md](algorithms/copilot-validation.md)):

- **Superstructure test calibrated.** `analyze_modulations` called a zone-boundary point ordered
  when it was in the top 1 % of the modes, with no correction over the 7 points × elements tested,
  and found "order" in 6 of 12 configurations without any (pure thermal motion included). A point
  now needs ratio ≥ 3 and ≥ 3σ against independent displacements (χ² with the sites' degrees of
  freedom, Šidák-corrected); each zone-boundary row carries its σ. Planted order stays > 8σ.
- **Cluster check no longer denies a distortion its own axis test shows.** "LESS distorted than
  independent atoms … no correlated shape distortion on top" was printed for Jahn–Teller octahedra
  whose <100> elongation the gyration-axis test found; the sentence now names the favoured axis.
- **Two answer checks.** "Polar" for a non-polar group (the local-symmetry skill promised it), and
  a superstructure claimed when every modulation analysis found none. A JSON key ("chiral": false)
  is no longer read as a claim.
- **The validation suite** (`web_app/frontend/validation/`): 9 kinds of planted-truth case, each a
  dev and a held-out test case; arms A (default pages), A+ (engine sweep), B/C (the Copilot with
  and without skills) and D (truth) scored by one fixed function; transcripts, number provenance,
  a spend cap, and the supervisor findings schema. Arm A: 0.81 on both sets; it reads ordered
  antiphase displacements as random off-centring, misses cluster breathing and rigid dimers, and
  the 0.2 Å card misses a polar shift below its tolerance. The review's larger findings (split shells cut
  at a deep trough, the cluster null under collective motion) are listed in
  [copilot-review.md](algorithms/copilot-review.md).

### Fixes from the 2026-10-12 review

An independent review of main (c26f8d4) found these; each fix has a regression test built from the
reviewer's repro.

- **Clusters.** A site at fraction 0 whose mean sat a hair below the cell edge wrapped to 0.9999,
  and its clusters paired atoms of neighbouring cells (a rigid V₂ dimer: z −32.5 → −0.6, broken
  edges 2.3 %): templates, bases and outlier positions now use the unwrapped site mean (both
  engines). Templates of one kind are measured with their own bonds (two Mo₃ chains listed 0–1–2
  and 1–0–2 gave a ladder width of 1.03 Å). The bonds of a cluster missing an atom are no longer
  counted as contacts between clusters, and on a tie the analysed kind leads the list of kinds and
  never appears among the other sizes.
- **Modulations.** The percentile is taken over the compared modes, one per ±q pair without the
  long waves, with a rank: a cell doubling at 55× the median under a stronger box-scale wave was
  reported as no order.
- **Bond angles.** Only an A–B–A triplet in one window is read against an ideal polyhedron; a band
  with no angles says "none", not "null°"; 15° and 30° engine bins are shared out by overlap (no
  comb in the 10° histograms).
- **Answer checks.** A hypothetical, conditional or question mention ("to check whether it has
  converged", "until χ² plateaus") is not a claim, and a negation counts only just before the
  claim ("With no remaining drift in χ², the run has converged" now fires). A lost inversion is
  judged against the parent group too, and a subgroup the checks do not know is never pinned on
  the parent.
- **The Workbench's questions.** A window is asked about once per question and the choice applied
  to later calls; a question that names the shell ("only the first Ta–Se peak", "the whole first
  shell", "both peaks") asks nothing; "2.2 – 2.7 Å" stays one option; a typed window goes to the
  window asked about.
- **Context guard.** Tool results count at 1.4 characters per token (gemma4 reads them at 1.64,
  not the 2.5 assumed), a model LM Studio loads for the request is checked before its reply is
  shown, and Ollama's native path sends the API key.
- **Copilot view.** A step's result is what the model read (fitted to its budget); answering a
  question card keeps the focus on the card, then the composer; "(100)" in a 3D view is the plane's
  normal, not the direction [100].

### Docs: a one-screen README

The README is now the name and live link, what the Workbench does, a quick start, the AI Copilot in
brief, its status, a linked documentation index and the license. Nothing was dropped: the
screenshots moved to [QuickStart.md](../QuickStart.md#what-it-looks-like), the signature equations
to [ALGORITHMS.md](ALGORITHMS.md#signature-equations), the license notes and citation request to
[REFERENCE.md](REFERENCE.md#license-and-citation), and the Copilot details and its real-model
results were already in the [Copilot guide](../web_app/frontend/src/llm/README.md). A stale link
to `SaveMenu.jsx` (now in `src/ui/`) is fixed in run-dashboard.md.

### AI Copilot: skills (drafts), and a refusal for windows that cut a shell

- **Four skills, drafted** (`.claude/skills/<name>/SKILL.md`, shared with Claude Code, in MATERIA's
  format): `rmc-fit-assessment`, `coordination-shells`, `local-symmetry`,
  `displacement-anisotropy`. Each opens with "Draft: awaiting the user's interview" and holds the
  order of steps, judgment calls marked "(to confirm with the user)", an acceptance checklist, the
  known failure modes with the spot check or review that found them, and small reference files.
  `docs/SKILL_INTERVIEWS.md` holds 7–9 questions per skill for turning each into the user's method.
- **`read_skill(name, reference?)`** (`llm/copilot/skills.js`): the system prompt lists only the
  skills' names and descriptions (2,150 characters); a skill (4,700–5,500 characters) is read on
  demand with its own result budget and is never compacted later in the question (`resultBudget`
  and `keep` in the tool contract). The host bundles the files (`src/copilotSkills.js`, the
  `skills` prop; Vite `server.fs.allow`, Docker copies them); a malformed skill is reported.
- **`compute_bond_angles` refuses a window that cuts the coordination shell** (fewer than 60 % of
  the central atoms sharing the most common count, or most with none): the angles are withheld and
  the counts that show the cut returned, unless `window_confirmed` — set when the user keeps their
  window in the Workbench's question. The 60 % threshold is an interview item.
- **Replays** (`src/__tests__/copilotSkillReplays.test.js`): one per skill on the Demo with the
  real engines and a scripted model; the right behaviour passes, the old wrong one is refused or
  flagged.

### AI Copilot: the Workbench asks before the model guesses

With a real model (gemma4:26b, 2026-10-11) the Copilot never called `ask_user`: on a request that
left a choice open it picked one. Where the run data show that the choice changes the result, the
Workbench now shows the question card itself, before the model is asked or before the analysis
runs, with the options it computed and a default (`llm/copilot/clarify.js`).

- **Before the first request:** a bond-angle or coordination request naming only an element
  (which bonded neighbours: one option per partner shell, nearest first); a named pair or triplet
  whose first shell is split or runs into the next, when the window changes the coordination
  number by ≥ max(0.5, 15 %) (the whole shell, the first peak, or a cut 0.2 Å either side); the
  displacement directions of an element with several non-cubic sites (one site, or all pooled:
  pooling can cancel a per-site preference); the space group with no tolerance named, when more
  than one group holds clearly above the noise floor. The answer joins the user's question for
  the model.
- **Before an analysis runs:** a window that cuts a first shell (`compute_bond_angles`,
  `analyze_clusters`, through the new `shellWindows(args)` in the tool contract) is put to the
  user — the whole shell (default) or the window as asked — and the analysis runs in the one
  chosen.
- **Conservative:** a question that settles the choice (an ideal polyhedron one option's CN
  matches, a window or tolerance in Å, a site number) asks nothing; the six spot-check questions
  ask nothing. The questions are `ask_user` steps (`app: true`), count toward the two per request,
  and the model's own `ask_user` still works.

### AI Copilot: answers checked against what the analyses found

- **Checks on the answer** (`llm/copilot/answerAudit.js`): tools return `facts` (not sent to the
  model) and the run context adds the watchdog's status and the ladder's groups; an answer that
  says "converged" while the history is still improving, calls a group chiral or centrosymmetric
  against its point group, says an inversion centre is lost from groups without one, or claims a
  preferred direction no analysis found at ≥ 3σ gets a **Check** line under it.
- **The Outcome chip is lowered:** with any finding, an answer the model called *achieved* or
  *partial* (or gave no verdict) shows **Not verified**, with what the model claimed; *not
  achieved* stays. A computed contradiction outranks the model's self-assessment. gemma4's
  correct answers of 2026-10-11 raise no check.
- **Point-group facts** (`llm/context/pointGroups.js`): `detect_symmetry` and the context's
  tolerance ladder state each group's inversion centre, Sohncke class and polarity (the spot check
  called P-42₁m chiral).
- **Angles read correctly:** `compute_bond_angles` adds the sin θ-corrected 10° histogram and, for
  a whole shell of about six or four neighbours, the split against an ideal octahedron or
  tetrahedron (the Demo's TaSe₆: 80 % cis at 89°, 20 % trans at 158.8°, which gemma4 had read from
  the raw modes as "highly distorted, 75–105°"). Its description no longer says to pick windows
  around the first g(r) peak.

### AI Copilot: it asks instead of guessing (`ask_user`)

The spot check's worst answer came from a guess: asked for "Se–Ta–Se angles against an ideal
octahedron", a model picked the first-peak window and cut the TaSe₆ shell in two. When a question
leaves open a choice that changes the science, the Copilot now asks.

- **`ask_user`** (`llm/copilot/askUser.js`): one short question, 2–4 options with a one-line meaning
  each, *Other…*, and the agent's default. The request is normalized (plain strings read as
  options, `"label — description"` split, repeats dropped, the default matched by label or number)
  and a malformed one is refused with a message the model can act on.
- **The card** (`llm/components/QuestionCard.jsx`), in the docked panel and the bottom sheet:
  option buttons (1–4, ↑/↓), *Other…*, *Skip*; it takes focus unless the user is typing, the log
  keeps its question in view, and a collapsed dock marks the rail. Answered, it keeps the question
  and the reply; the turn footer does not count it as an analysis.
- **The loop waits** (`interactive: true` in the tool contract, `askUser` in the tools' context):
  never cached, no time limit, the wait off the turn's 10-min clock, at most two questions per
  request and one per step. The answer returns as the tool result ("use exactly this choice";
  a skip says "use your default and say you assumed it"). Stop, a run change or an unmount
  cancels it cleanly.
- **Prompt rules** (`ASK_RULES`): ask only when the answer changes which analysis runs or its
  parameters (shell or r-window, site or pair, tolerance, run) and neither the context nor an
  analysis can settle it; settle what the data can first; never on a clear question; offer a
  default; use the chosen option and say which.
- **Argument validation:** a bare string where an object with one required string field belongs
  is read as that field (`options: ["a", "b"]`).
- **With a real model** (gemma4:26b, 16 questions on the Demo, 2026-10-11): it never called
  `ask_user`, including on "What are the Se–Ta–Se bond angles?" — the run context now holds the
  whole first shell (2.2–3.23 Å, CN 6), which it used — and on open requests ("the bond angles
  around Ta") it chose the main polyhedron and said so. Clear questions took 6–37 s, none asked.
- **Tests:** replays on the Demo's split Ta–Se shell (an ambiguous question asks and the chosen
  window — option, default, typed, skipped — reaches `compute_bond_angles`; a clear question does
  not ask; cancel; limits; the turn clock), and the card end to end against a mocked server.

### AI Copilot: a prompt is never cut to fit the model's context

A local server runs a model with a fixed context window and shortens a longer prompt without an
error, so the model answers without part of its rules or the run context. The Copilot's first
request on the Demo is ~28,400 characters (~8,600 gemma4 tokens), and Ollama's OpenAI-compatible
endpoint cannot carry a context length: it uses the server default, which Ollama 0.40 picks from the
GPU memory it finds — 262,144 on the maintainer's Mac (nothing was cut: no `truncating input
prompt` line in its logs since 2026-08-26), but 4,096 on a start where it saw no GPU (the same Mac,
2026-09-27). Checked with gemma4:26b on the Demo's first request (8,437 tokens): at a 4,096 window
Ollama answers 200, drops the run context, cuts the rest to 2,051 tokens, and the model calls the
wrong tool; at 10,240 or more it reads everything.

- **`provider/contextGuard.js`** (`streamCopilotChat`, used for every Copilot request): for an
  Ollama server it sends the request through the native `/api/chat` (`provider/ollama.js`, its own
  NDJSON parser) with `num_ctx` = the loaded context when it holds the prompt and 8,192 tokens of
  reply, else the next multiple of 8,192 (≥ 16,384, kept for the conversation, ≤ the model's
  maximum), and `num_predict` inside it, so a long reply cannot push the prompt's start out of the
  window. The prompt size is a deliberately high estimate (2.5 characters per token; ~32 % above
  gemma4's count on the Demo).
- **Refuses rather than cuts:** when even the model's maximum cannot hold the request, or LM
  Studio's loaded context cannot, the turn stops before anything is sent and says what to change.
  Cloud providers and other servers are unchanged; only `localhost` or port 11434 is probed.
- **A backstop:** Ollama's `prompt_eval_count` is what the model read; under 0.4 of the estimate
  (an uncut request reads 0.74 of it) the turn stops with an error.
- **A reply cut at the length limit** says so in the turn's footer.
- The setup guide, README and QuickStart no longer ask for `OLLAMA_CONTEXT_LENGTH`.

### UI kit: no sideways scroll from closed info popovers

At iPad-portrait widths (820–834 px) the Displacement Directions page scrolled sideways: a `?`
badge in the wrapped controls bar sat within 20 rem of the right edge, and its closed popover,
hidden with `visibility` only, still counted toward the page's overflow. A closed `InfoBadge`
popover now takes no space (`display: none`, still fading in where display transitions are
supported), and an open one takes the other edge when its own would run off the window
(`ui/popoverEdge.js`). Checked on every page at the kit's device sizes.

## 0.7.0 — 2026-10-11

### AI Copilot: docked beside the pages, pointing at what it explains

The Copilot no longer floats over the pages, and it now shows what it explains on the page itself,
with the design language of the maintainer's other agents (MATERIA, NEBULA3D), re-implemented here.

- **Docked, not floating.** On a window ≥ 1080 × 541 px (laptops, monitors, tablets in landscape)
  the Copilot is a side panel in its own
  grid column (`App.css`: `.app-container` gains a column and a row; `.app-copilot`), so it never
  covers a card. The pages shrink and reflow beside it: the header's and the pages' width tiers
  (1599/1040 header rows, the 980/1100/1200 page stacks, the 1500 three-column dashboard, the
  footer, `.ui-kbd`) are now **container queries** on `.main-content` (`container: workbench`),
  so a 1440 window with the Copilot open lays out like a 1040 one. The phone tiers (≤ 760, ≤ 540
  tall) stay media queries. Drag the left edge to resize (keyboard: ←/→, Home/End; double-click or
  Enter resets); the pages keep ≥ 760 px. On tablets in portrait, phones and short windows it is a bottom sheet
  in its own grid row (resizable from its top edge; `--copilot-sheet-h` is taken off
  `--viewport-h` inside the pages); when the agent brings a card on screen the sheet comes down,
  as far as 28 % of the window, until the card fits whole above it. Collapsed, it is a slim rail (a bottom bar on a sheet) with a
  red dot when an answer lands while collapsed. Width, sheet height and open state persist
  (`rmc-llm-dock-v1`). The floating launcher, the *Larger* toggle and the dialog role are gone;
  the dock is an `aside` landmark.
- **Pointing at cards** (`highlight_cards`, a new view tool built from the host's card catalog,
  `src/copilotCards.js`). Each Workbench card carries `data-copilot-card` (the model card, Detected
  SG, the χ² history, one per plot file, and every card of the analysis pages). The agent's cards
  breathe a blue edge (2.4 s; a steady edge under `prefers-reduced-motion`); the main one's page
  opens and it scrolls into view. On a plot it can mark x positions — the g(r) peaks it names —
  drawn as dashed blue guides with ▽ and a label, and the plot zooms to them (`InteractivePlot`
  `marks`). The card of a running analysis (`assess_convergence`, `find_gr_peaks`,
  `compare_fit`, …) breathes softly. A highlight ends at the next one, a new question, a click on
  the page, or when the Copilot collapses.
- **Card chips.** Answers cite cards as `[label](card:<id>)` or `[[id]]`; each renders as a chip
  that brings the card back (from another page too). Unknown ids render as plain text; the
  `card:` scheme is let through react-markdown and the sanitizer for chips only.
- **3D views** (`orient_3d_view`, `current_view`; `src/copilotViewTools.js`). The six 3D views
  (Atomic Density's folded cell, Bond Geometry's bonded cell, the PCA density and site ellipsoids,
  the displacement sphere and its site picker) join a registry while mounted
  (`src/viewer3d.js`, `src/useCopilotViewer.js`). The agent can look down `[u v w]` (any
  spelling: `[111]`, `[1-10]`, `[1 -1 0]`, `[0.71 0.71 1]`) or a, b, c, PC1–PC3, with c up (b up
  looking down c), in oblique cells too; centre a site, ring sites in breathing blue, select a site
  as a click would (waiting for its map or volume to load), or reset. A turn is a 650 ms slerp
  that lands exactly on the pose (`src/cameraPose.js`). On the sphere's PCA frame a crystal
  direction is projected onto the site's principal axes. `current_view` reports the page, its
  cards, the highlight and each 3D view's direction (as the simplest `[u v w]` within 2°) and site.
- **The rule changed:** analyses still only compute; the **view tools** (`view: true` in the tool
  contract, never answered from the duplicate-call cache) change only what is on screen — never
  data, an analysis input, the shared tolerance or a file. The prompt (`COPILOT_PROMPT_VIEW`) asks
  the agent to point before it explains, cite cards as chips, and turn the 3D view to the direction
  it describes.
- **Tool contract.** `run(args, { signal, focus })` (the dock's card-focus controller),
  `cards(args)` (cards to breathe while running), `view`; argument schemas gain arrays and objects
  (`items`, `minItems`, `maxItems`, nested `properties` / `required`), and a list sent as JSON text
  or a comma list is accepted. `CopilotChat` gains `cards`, `onMarks`, `onLayoutChange`,
  `className`, `storageKey`.
- **Checked** in Chrome with Playwright at 1440 × 900, 1920 × 1080, iPad Pro 13" landscape
  (1366 × 1024) and iPad landscape (1180 × 820), docked, and iPad portrait (820 × 1180) and a phone
  (390 × 844), as a sheet: open and collapsed, every page with the Copilot open, a highlighted card in view
  and not under the Copilot, a chip, a click that clears, a [111] view, and the remembered state.
- **Checked with a real model** (gemma4:26b, Ollama, the Demo, 2026-10-11): five questions —
  a highlighted g(r) card with a mark at the Ta–Se peak, a chip that opened PCA Ellipsoid, the
  direction sphere turned down [111], site 1 centred and ringed, and `current_view` — all
  passed with no page error. Two display bugs it found are fixed: the Outcome note now renders
  its math, and a chip written as the bare reference shows the card's title
  (`src/llm/README.md`, *The view tools with a real model*).


## 0.6.0 — 2026-10-10

> This work was first tagged v1.0.0 on 2026-09-26 and the tag was withdrawn the same day; the
> work continues toward 0.6.0, which will also bring the Bond Geometry redesign and device layouts.

Every engine's math and physics was audited end to end, and the audit's roughly 100 defects were
fixed in both runtimes. They included wrong numbers, Python ↔ browser disagreements, inputs that
crashed or silently produced garbage, and docs that described code that no longer existed. Nine
groups fixed them test-first (Auto StoG ×2, structure KDE, PCA ellipsoid, displacement directions,
bond angles, parsers and dashboard, symmetry finder, Flask API). Each group was reviewed
independently, then merged and integrated. Several numbers the app reports change *because they
were wrong*: read **Upgrading from 0.5.0** before comparing 0.6.0 output with 0.5.0.

### AI Copilot replaces the AI Assistant

The **AI Copilot** replaces the chat-only AI Assistant. It is a **chatbox over every page**, not a
tab: an *AI Copilot* launcher in the bottom-right corner opens a floating panel. The panel can be
made larger and closes with ✕ or Esc, keeping the conversation. On a phone it is a full-height
sheet. The Assistant could only answer from one fixed summary of the run. The Copilot also **runs the
Workbench analyses itself** through OpenAI-compatible tool calling, reads what came back, and says
whether that answers the question. The plan and its review are in
[COPILOT_PLAN.md](COPILOT_PLAN.md); the method is in
[ai-assistant.md Step 17](algorithms/ai-assistant.md#step-17--the-copilot-tool-loop-analyses-checks-budgets-verdict).

- **Ten analyses.**
  - Over the run data: `read_run_data`, `assess_convergence`, `find_gr_peaks` (peaks and the first
    coordination shell), `compare_fit` (where the model misses the data).
  - Through the same engines the pages use, in the browser or through Flask: `list_sites`,
    `analyze_displacement_directions` (calibrated significances, the peak direction as `[u v w]`),
    `compute_bond_angles`, `analyze_clusters` (cluster distortion against independent atoms, a new
    engine), `analyze_modulations` (order that enlarges the cell, a new engine) and `detect_symmetry`
    (any tolerance, without changing the Detected SG card).

  A tool is offered only when the run has its inputs. The tools only compute: none switches a
  page, changes the shared tolerance or writes a file.
- **Show in \<page\>.** A step's button switches the page *behind* the chat and loads the same
  analysis there:
  - Bond Geometry fills the triplet and windows and computes; on the demo it matches the
    Copilot's 41,336 Ta–Ta–Ta angles and coordination 2.73.
  - Displacement Directions selects the site and the resolution ('auto' when the Copilot let the
    engine choose). The page keeps its own smoothing, so the drawn map can differ from the
    Copilot's unsmoothed one; the significance tests use raw counts either way.
  - Other steps just open their page.

  The button is the user's action, through an App-level page-request channel
  (`pageRequest` → the pages' `request` prop, applied once per click; `tripletInputsFromRequest`
  in `workers/triplets.js`). On a phone the sheet closes so the page is visible.
- **Checked results.** Every result carries deterministic checks computed from the full engine
  output, e.g. "no angles in these windows — widen them", "no significant peak direction
  (1.2σ < 3σ)", or "the first g(r) peak sits on the closest-approach constraint". A failed check
  goes back to the model, which changes its arguments and retries.
- **A verdict.** Every answer ends with `Outcome: achieved | partial | not achieved — <why>`, shown
  as a chip. It is the model's self-assessment; the checks on the step cards are the computed part.
- **Step cards.** Each card shows the status, the duration and the most severe check. Expanded, it
  shows every check, the arguments and the exact data the model read. A polite live region
  announces step and verdict changes.
- **Bounded.**
  - Per question: ≤ 6 rounds of tool calls, ≤ 4 calls per round, 120 s per call, no new round
    after 10 min.
  - Each result is fitted to 3,000 characters and a question's results to 12,000, so a small local
    context window keeps the run context. Older results are compacted to their checks first.
  - An identical repeated call returns the cached result.
  - When the budget is used, a last request with `tool_choice: 'none'` asks for the answer.
- **Robust across providers.**
  - Streamed tool calls are assembled by index, id or order.
  - Gemini's `thought_signature` (`extra_content`) is echoed back, and whitespace-only assistant
    content is sent as `null` (Anthropic).
  - A model that "does not support tools" (Ollama) falls back to plain answers for the session.
  - An empty last round (Ollama ignores `tool_choice`) closes with a fixed answer.
- **Stop cancels the work.** A browser-loaded run is analysed in a worker of the Copilot's own,
  terminated on Stop or timeout, so the pages' shared worker never queues behind it. Server-folder
  requests are cancelled.
- **One run per conversation.** Opening another run restarts the conversation. Earlier questions
  travel as question + final answer only, and an unanswered question is dropped.
- **Settings.** The drawer adds *Let the Copilot run analyses* (default on). The cloud warning now
  says that analysis results are sent too. The setup guide asks for a tool-capable model and a
  larger Ollama context (`OLLAMA_CONTEXT_LENGTH=16384`).
- **Code.**
  - New in `src/llm/`: `copilot/` (agent loop, tool-call assembly and validation, tool contract,
    run-data tools, Outcome parser, turn reducer), `prompts/copilot.js`, `useCopilot.js`,
    `CopilotChat` / `CopilotView` / `ToolStep` / `MarkdownMessage`.
  - The host builds the engine tools in `src/copilotTools.js` and passes them as a prop, so the
    module still imports nothing from the app.
  - `useSiteCloud.js` exports its routing as `requestAnalysis()`, with an optional worker and abort
    signal. Request params can no longer override `kind`, `text` or `dir`.
  - `provider/client.js` sends `tools` / `tool_choice` only when given, yields tool-call fragments,
    and attaches the HTTP status to errors.
  - Removed: `ChatView`, `AssistantPage`, the unrendered `AssistantPanel` / `AssistantWorkspace`
    card, `useStreamedReply`, and `buildChatMessages`.
  - The Copilot tab and its page are gone. The Dashboard's Copilot-only work (symmetry ladder,
    `.dat` parse) starts when the chatbox first opens.
- **Tests.** New and extended vitest cases cover:
  - the agent loop against a scripted model;
  - SSE tool-call parsing, including the Gemini shape;
  - argument validation and the Outcome parser;
  - the run-data tools;
  - the host tools on the committed demo run with the real engines (F-43m at 0.2 Å; the Ta₄
    cluster's Ta–Ta–Ta mode near 60°; an empty window fails its check);
  - the chatbox end to end against a mocked server (nothing fetched before it opens, Show in
    hands over the page parameters, Esc keeps the conversation, a run change restarts it);
  - `tripletInputsFromRequest` round trips.
- **Tested with real local models** (2026-10-10, gemma4:26b on Ollama, qwen3.8-27b on LM Studio,
  six questions on the Demo; results in
  [`src/llm/README.md`](../web_app/frontend/src/llm/README.md#tested-with-real-models)). Tool
  calling worked with no malformed call. The runs found four bugs, fixed:
  - a verdict line in backticks (as the prompt shows it) or under a heading marker now gives its
    chip, instead of staying in the answer text;
  - `find_gr_peaks` notes when the second peak overlaps the first-peak window (on the demo, Ta–Se
    at 2.52 and 2.74 Å: the old window cut the TaSe₆ shell) — since replaced by the first-shell
    window below;
  - `compute_bond_angles` warns when fewer than 60 % of the central atoms share the most common
    bond count (`MODAL_COORDINATION_SHARE`): the window may cut through a shell;
  - an error event sent mid-stream (LM Studio out of GPU memory) is shown, instead of ending the
    turn empty and crashing the chat (`turn.items is not iterable`); an empty reply shows its note.
- **A scientific review with a real model** (2026-10-10, gemma4:26b, the Demo and the 5 K GaNb₄Se₈
  run). The answers were fluent and the tool calls valid, but several conclusions were wrong for
  reasons the tools and the context could fix:
  - *"Converged — a plateau."* The 48-point history hid a still-falling tail. The context now
    carries the watchdog's own verdict (`watchdog_status`: improving) and the change over the last
    tenth of the run (−0.11 in ln χ²).
  - *"Two datasets."* The F(Q) and its Fourier transform are one measurement. The context says so
    (`fitted_in`, `independent_datasets`, `same_measurement_as`) and adds the **Faber–Ziman weight of
    every partial** per radiation: under X-rays Ga carries ~5 % of the GaTa₄Se₈ signal.
  - *"No first peak is constraint-limited."* The Ga–Se distribution is cut at the 2.2 Å constraint
    (g there is 9 % of the peak), which the half-maximum test missed. `find_gr_peaks` now checks g at
    the constraint and analyses the **first coordination shell**: its window to the first deep
    minimum (split shells kept whole), mean distance and spread from r²g(r), coordination numbers
    both ways with their sensitivity to the cut, and whether the shell is separated from the next
    (demo Ta–Ta: CN 2.76 of 3, distances filling the 3.3–4.0 Å gap). The context's pair
    correlations carry the shell's mean and CN too.
  - *"The Ta–Ta peak is shifted to 2.93 Å: the cluster is distorted."* A g(r) peak sits below the
    mean of a broad shell (3.00 Å), and its width cannot tell a distortion from thermal spread. The
    new **`analyze_clusters`** tool ([ai-assistant.md Step 18](algorithms/ai-assistant.md#step-18--cluster-distortion-against-a-shuffle-null-analyze_clusters);
    `rmc_toolkits/clusters.py`, `workers/clusters.js`, `/api/clusters`, a parity fixture) measures
    every cluster against a shuffle null of independent atoms: the demo's Ta₄ clusters are 23 %
    *less* distorted than independent atoms (rigid units), favour compression along <111> and,
    weakly, elongation along <110>. The checks require an effect size as well as a z-score.
  - *"No preferred direction, so no cluster distortion."* The model pooled all 16 Ta sites, whose
    local axes differ. A pooled map now says so and points to single sites or `analyze_clusters`.
  - *"Ga–Se 2.428 ± 0.02 Å."* The ± was invented. `compute_bond_angles` now returns the bond-length
    sd, and the prompt forbids an uncertainty no tool gave.
  - *"Where does the fit fail?" — no tool.* The new **`compare_fit`** tool gives Rw and the share of
    the misfit per range of r or Q and the largest residual features (demo G(r): 62 % of the misfit
    at 2.3–5.3 Å).
  - The symmetry ladder's lower rungs were read as distortions. The context now gives the ladder's
    **noise floor** (3√2 · max rms displacement / √copies; 0.031 Å on the demo, where F-43m first
    holds) and says the folded average cannot show a superstructure.
  - *"Could a superstructure be hidden?" — no tool could say.* The new **`analyze_modulations`**
    ([ai-assistant.md Step 19](algorithms/ai-assistant.md#step-19--displacement-modulations-order-that-enlarges-the-cell-analyze_modulations);
    `rmc_toolkits/modulations.py`, `workers/modulations.js`, `/api/modulations`) Fourier-transforms
    every site's displacements over the cell grid and tests the zone-boundary points against the
    configuration's own diffuse background. Neither real run holds a cell doubling, though GaNb₄Se₈
    is orthorhombic with a doubled cell at 5 K; both carry box-scale displacement waves.
  - The context is serialized compactly (arrays and small objects on one line) under a 6,000-
    character budget, so both real runs fit untrimmed.

### Detected SG: download the symmetry-averaged structure as a CIF

A **`⤓ CIF · <group>` button** under the Detected SG title downloads the selected group's structure
as a CIF. Pick a group by clicking its brick on the tolerance ladder, as before, and the button
follows. The button is placed like Model information's file name, so the card keeps its height. The
CIF is the RMC model folded into one unit cell, averaged in the group that brick stands for, and
written in the standard cell the symbol is named in
(`averageStructure.js`, `cifWriter.js`, `symmetryCif.js`; the method is in
[run-dashboard.md Part C](algorithms/run-dashboard.md#part-c--the-symmetry-averaged-cif-of-a-ladder-brick)).

- **Averaged, not idealized.** Each reference site's arithmetic mean over its box copies is
  averaged over its orbit through the group's operations, so special positions are exact and free
  coordinates keep their measured values (GaTa₄Se₈ demo in F-43m: Ta and Se on 16e at x = 0.6463,
  0.3854, 0.8856; Ga on 4c).
- **Exact operations.** The finder's least-squares translations are made into an exact group by
  averaging the product 2-cocycle. The origin then moves by the smallest shift that puts every
  operation on the 1/48 grid (0.0005 Å on the demo; a box cut on an arbitrary origin gets an origin
  on its symmetry elements). The shift, the cell transformation and the largest symmetrization shift
  are recorded in the file's header, and the operations are listed explicitly.
- **U_ij from the atom distributions**: the pooled second moment of every atom of the orbit about
  its symmetrized position, in the site-symmetry form. A distortion the picked group averages away
  shows up here (pick F-43m on an R3m-distorted model and the Ga U grows by the Ga offset).
  Mixed sites get per-element occupancies.
- The cell metric is averaged over the point group. A group with no standard cell (a crystal class,
  a `≥` lower bound) is written in the `.rmc6f` cell with no symbol. Orbits are taken at the
  brick's tightest tolerance, not at the midpoint the click selects.
- `structureFromRmc6f` basis sites now also carry `count`, `elementCounts`, `mean` (arithmetic,
  unwrapped about the circular mean) and `covFrac`; the finder still reads `frac`.
- A failed export is listed in the page's Problems section (`CIF · <group>`) until an export
  succeeds.
- **On ITA's origin.** For a named group the structure is moved onto International Tables' own
  origin, the equivalent one nearest the `.rmc6f` origin (never along a polar axis). The shift is
  solved exactly from the integer congruences (`itaOrigin.js`). The CIF then lists ITA's operations
  and spells origin choice 2 (`F d -3 m :2`), so the H–M symbol alone describes the file, and every
  site gets its Wyckoff letter. The demo's Cmm2 and P-4̄2₁m exports had several `?`; now none do.
  The demo's F-4̄3m box already sits on an ITA origin and does not move (Ga on 4c); diamond built on
  origin choice 1 lands on choice 2's 8a. The operations come from a new table, `itaOperations.js`
  (230 generator sets, the settings of the Wyckoff table). The table equals spglib's operation sets
  exactly (`tests/generate_ita_operations_fixture.py`, spglib dev-only), and the symmetry tests'
  fixture now reads it, so there is one copy. A group with no standard cell keeps the 1/48-grid origin
  and says its listed operations are authoritative.
- Flask server-folder mode still has no Detected SG card.

### UI kit (no visual change)

Every workspace page now draws its chrome from one shared UI kit, `web_app/frontend/src/ui`
(`ui.css` plus thin React components; see its README), instead of per-page look rules — the
polished look had lived in `PcaKdePage.css` as `pca-*` classes that other pages borrowed by name.
Theme-invariant values that were hard-coded (brand blue, status and plot-kind colours, chart ink,
radii, durations, control heights, the type scale) are tokens in `index.css`. Page stylesheets
keep layout only, apart from a few listed exceptions (the app header and brand lockup, the PCA
statistics dividers; see the kit README). `InfoBadge` and `SaveMenu` moved into `src/ui`; the footer is
`components/AppFooter.jsx`. The refactor was verified pixel-identical against `main` on every page
(desktop, laptop and phone; light and dark) — the first viewport and, because each page scrolls
inside its own root rather than the document, every below-the-fold scroll step — with computed
styles, hover and keyboard-focus states included, compared page by page.

### Lean pages

Page bodies no longer carry explanatory sentences. Explanations and caveats moved, intact, into the
nearest `?` help (the kit's `InfoBadge`) or a title tooltip; state flags became compact cues; only
what a user needs to act stays on the page — one-line errors, validation messages, the no-run
prompt and one short empty-state line per card. The same kind of text gets the same treatment on
every page:

- **No run:** one line everywhere — *"Open a run folder with an `.rmc6f` file."* on Atomic
  Density, Bond Geometry, PCA Ellipsoid and Displacement Directions (Atomic Density showed a red
  error), *"Open a run folder."* on the Dashboard and the AI Assistant.
- **Header:** the startup "Files stay on your device" notice is a *Where your files go* `?` beside
  the run-folder field; the non-Chromium Live Data notice is one line.
- **Dashboard:** the status banners are gone (a first-load *Loading…* covers Flask); a Live Data
  read of a half-written `.rmc6f` shows a *previous read* chip in the Model information heading;
  other runs in the folder are a *+N runs* chip; several failed χ² logs are one counted line with
  the per-log list behind *Error details*; the parse warning's sentence is in the card's `?`.
- **Atomic Density:** the notes under the KDE canvas are one short line per reason, keyed on the
  new `messageCode` (payload sentence in the tooltip), and the reasons — declined slabs, grid
  resolution, the 3 : 1 anisotropy caveat — are in the KDE Slice `?`. The header shows `browser` or
  `server`; the help no longer claims the browser on the server path. The canvas no longer prints
  *"No atoms in this slab"* before the first result or after an error.
- **Bond Geometry:** the two hints fold into the angle card's placeholder (*"Pick a triplet, then
  Compute."* / *"New configuration — Compute again."*) and the triplet / window helps; the
  Partial PDF placeholder names the missing pair when the partials file lacks it.
- **Auto StoG** (hidden in the shipped build): the *How Auto StoG works* explainer, the intro
  and the fieldset descriptions are `?` helps on the cluster labels and legends; the readout cards
  have short subs with their reasons in label helps (the Concordance card turns neutral when
  a_fz is unresolved); long engine refusals show a one-line summary with the full
  message behind *Error details*.
- **PCA Ellipsoid / Displacement Directions:** *"Atoms skipped: …"* banners and a short grey-shell
  legend item (full sentences in their tooltips); a zero-spread site reads *No displacement.*, its
  reason in the principal-axes `?`.
- **AI Assistant:** drawer and chat hints are tooltips or part of the *How to connect a model*
  guide, which a failed Test now opens; the unreachable-server hint is one line pointing to it.

Additive data fields back the short forms; no message string changed: `messageCode` on every KDE
slice payload (both runtimes, `/api/kde/slice` included), `summary` on the browser Auto StoG
engine's diagnostic errors (and the worker's error reply), `parseErrors` beside `parseError` on the
static-mode χ² history, and `unreachable` on a failed AI Assistant connection check. The kit's
`ui-intro`, `ui-disclosure`, `ui-fieldset__desc` and `ui-controls__note` lost their last user and
were removed.

### Problems section on every page

- **Errors are never silent.** Every page lists what is wrong on it in one compact section at its
  top, errors first, then warnings: a site table or structure that failed to load, a KDE or Compute
  that failed, atoms the parser skipped, a plot file or χ² log that failed to parse, a figure that
  could not be saved, and any error nothing else caught. Collapsed it is at most two lines (two
  problems, or the first and "N more"); a repeated failure counts up (×N). Nothing appears when a page has no problem, so the
  layouts are unchanged (checked pixel for pixel on every page with the demo run).
- It replaces the page-level danger and caution banners (PCA Ellipsoid, Displacement Directions,
  Bond Geometry, Atomic Density, the Dashboard's run-folder and structure errors, the save-failure
  banners). A card keeps its own cue where it explains an empty view (a failed chart, a 3D error
  badge, the Bond Geometry prompt).
- **Newly surfaced:** failed saves of the 3D views and canvases (PCA view, Site ellipsoids, the
  direction sphere, KDE slice, slab, folded cell, Bond Geometry's bonds view), which used to fail
  without a word, and a capture that returned no image; Bond Geometry's skipped-atoms warning, which
  that page never showed.
- A failed save stays listed until the same save succeeds, it is dismissed, or another run opens.
- **Fixed:** in the static build with no run open, PCA Ellipsoid, Displacement Directions and Bond
  Geometry showed "Request failed with status code 404" (they asked the server that the static
  build does not have). They now ask nothing until a run is open.

### Atom colours: one map per model, similar colours separated

- **One element-colour map on every page.** Atomic Density, Bond Geometry, PCA Ellipsoid and
  Displacement Directions all build it from the model's full species list and atom counts, so an
  element is drawn in the same colour everywhere (checked on the demo run).
- **Similar colours are separated automatically.** Two table colours can look alike (Ti/V, Fe/Co,
  Fe/O, Mn/Sn, Cr/Mn). When two elements in a model are closer than 0.15 in OKLab, or than 0.06 for a
  simulated deuteranope or protanope, the less abundant one moves to the nearest colour that is far
  enough from all the others. The majority species keeps its familiar colour (Fe₂O₃ keeps O red and
  moves Fe). A model whose colours are already distinct is unchanged: GaTa₄Se₈, GaNb₄Se₈ and BaTiO₃
  keep every table colour. Hovering a moved element in a legend says what it was too close to.
- **Ta has its own colour** (`#F39B7F`, NPG salmon) instead of the fallback cyan next to Se's teal.
- Details: [structure.md, Step 8](algorithms/structure.md#step-8--element-colours-shared-by-the-slab-canvas-the-3d-view-and-the-legend).

### Bond Geometry Phase 1

The first phase of the Bond Geometry redesign: presentation only, no engine change (the payload,
the CLI and every number are as before). Details in
[bond-geometry.md](algorithms/bond-geometry.md#bond-geometry--the-page).

- **Hero layout.** The angle distribution is the hero card (7 of 12 columns, full height); the
  folded cell and the partial g(r) stack on the right; one column below 1100 px. The grid takes
  the height left under the model card and the controls, and the bar keeps one row down to 1440 px
  even with the split window on, so there is no page scroll at 1440 × 900, 1600 × 900 or
  1920 × 1080. The *Triplet result* strip became a KPI rail inside the hero — angles per
  central atom, mean coordination with the modal share, mean bond length per bond type — which
  reads "—" before Compute, so nothing moves when a result lands. The mean ± std angle moved off
  the headline to the Angles tile's sub line (a multimodal mean is not a bond angle).
- **An angle axis.** 0–180° unpadded, labelled every 30°, minor marks every 10°, vertical grid at
  the labels, y from 0; the histogram draws as a step curve with a light area. No reference line:
  a dashed *random bonds* line (1 when sin-corrected) was tried and dropped, as the axis label
  already says random = 1 and against a crystal's peaks (6–22× on the demo) it read as a stray
  gridline. The axis labels say what they are
  (`angle at Ta, θ (°)`, `sin-corrected (random = 1)`). These are opt-in `plotData` fields of
  `InteractivePlot` (`xDomain`, `xTicks`, `xMinorStep`, `xGrid`, `yMin`; per series `curve:
  'step'`, `fill`, `width`, `legend: false`); every other plot renders unchanged
  (`interactivePlotAxes.test.jsx` pins them to the markup from before).
- **One colour system.** Element chips (a dot, a tint, a ring on the central atom) name the
  triplet in every card title; `BOND_COLORS` (`plotPalette.js`: A–B = plot blue, B–C = plot
  orange) colour the 3D bond sticks (previously two unrelated hard-coded colours), the split
  window guides and chips, the partial curves and the bond dashes; references stay neutral grey.
- **States inside the cards.** Before Compute the hero shows the axis dimmed with a prompt card
  (triplet, window, Compute, the central-atom count, the box, where it runs); computing shows a
  shimmer or dims the shown plot, with a centred badge and a spinner on the button; a result whose
  inputs changed shows an *inputs changed* chip and the button reads **Update**; a validation
  error is the prompt's line and marks and focuses the field it names; a result with no triplet
  in the windows keeps the empty axis with "No A–B–C triplets in these windows." instead of a flat
  zero curve; a Live Data configuration change shows a *new configuration* chip. Without `PDFpartials.csv` the g(r) card is a slim row.
- **Controls as a form.** Enter in any field computes (an *Enter* key hint beside the button from
  1500 px up); Compute sits at the right with a ▶ glyph; numbers carry their unit inside the field
  (`[2.00 Å]`, `[1.0 °]`); the window labels are the bond roles (*A–B*, *B–C*) with their
  bond-role bar; a ⇄ button swaps A and C.
- **Titles and legend.** "(Ta)–(Se) partial g(r)" with a live window chip and the plot's Save in
  its header — the title names the curves, so the g(r) plot has no legend row and keeps that
  height (InteractivePlot's opt-in `legend={false}` and `actionsTarget`); the folded cell and the
  g(r) share the right column equally. "(Ta)–(Se) bonds" for the 3D card, whose element legend
  moved into the canvas as a pill with a bond swatch per window.
- **UI kit.** New pieces, documented in `src/ui/README.md`: `ElementChip`, `BondDash`, `KpiRail`,
  `Kpi`, `UnitField`, `PrimaryButton run` (`busy`, `stale`), `ControlsBar as="form"`,
  `InfoBadge side="above"`, and the classes `ui-legend--overlay`, `ui-overlay-center`, `ui-dim`,
  `ui-shimmer`, `ui-prompt`, `ui-role-bar`, `ui-kbd`, `ui-reveal`; all motion is off under
  `prefers-reduced-motion`.

### UI fixes from the lean-pages sweep

- **Failed figure saves say so.** A chart's *Save* failure (e.g. *Could not rasterize the figure*)
  never showed, and the Dashboard's *Save all figures* dropped its failures unhandled. Both are now
  a one-line danger banner under the chart toolbar / the Loaded-files header, cleared by the next
  save (the Save-all one also by opening another run folder).
- **Dashboard, server unreachable:** a failed run-folder listing with no server message reads
  *Could not list the run folder* instead of leaving only *Open a run folder.*
- **Displacement Directions:** with no site, the sphere's stage reads *Loading sites…*, then *No
  site selected.* (as the PCA Ellipsoid table does) instead of a blank canvas; the *About site clustering* `?` now explains the
  count/copies figure like PCA's.
- **PCA Ellipsoid:** the Site-ellipsoids panel shows its *Loading…* chip while the sites load.
- The never-shown *Structure data loads when needed* status is gone from static-mode runs
  (`structureError` is `null` when a `.rmc6f` is found).

### Device layouts

The workbench is laid out for the screens it is used on: 1080p ("1K"), 1440p ("2K") and 2160p
("4K") monitors, current Macs, iPads and iPhones. Every page was checked in headless Chromium at
the browser content size of each target (window chrome subtracted), before and after.

- **2K / 4K monitors: the UI scales.** The root type size — and with it every rem-sized control,
  gap, card and height budget — is 15px up to a 1080p-class window and grows linearly past it:
  17px at 2560 × 1440, 21px at 3840 × 2160 at 100% OS scaling, following the smaller of the
  width and height terms (an ultrawide or portrait monitor gets the size its short side supports).
  On a 27" 1440p or a 32" 4K panel, text keeps roughly the physical size it has on a 24" 1080p
  one instead of shrinking with the pixel pitch; 4K at 150% scaling reads as 1440p. Laptops,
  tablets, phones and 1080p monitors keep 15px (at 1920 × 960 every page renders as before; the
  one visible difference is the centred "Computing KDE…" message below). What drew in CSS pixels follows the same
  scale (`uiScale.js`): the Atomic Density canvas labels, and the Bond Geometry plots, whose
  viewBox now counts one unit per `uiScale()` pixels. Fixed pixel sizes that would have stayed
  small became rem (dashboard plot height caps, the loaded-files grid, range sliders, the
  tolerance ladder, the AI Assistant column). On 4K the Displacement Directions and Bond Geometry
  rows are no longer stopped short of the footer by their rem caps (now `max(60rem, 44vw)` /
  `max(56rem, 38vw)`).
- **MacBooks and 1366-wide laptops: one header row.** The full header needs ~1535px, so at every
  current MacBook's default scaling (Air 13" 1470, Pro 14" 1512) and on an iPad Pro 13" in
  landscape it wrapped to two rows. Between 1041 and 1599px it is tightened (tab and cluster
  spacing, the field bar's micro-label visually hidden) to ~1365px and stays on one 64px row
  (the Flask build's wider run-folder field fits from ~1410px).
  The Model information / Triplet result rails no longer drop their stats under the title cell
  when they wrap at laptop widths: the stats wrap beside it. The PCA Ellipsoid page no longer
  draws its 3D card over the footer on a 13-14" MacBook: side by side on a screen under 960px
  tall the canvas floor is 16rem (22rem elsewhere), the layout never shrinks below its content,
  and the WebGL canvases sit out of flow so their last pixel size cannot hold a card open.
- **iPad (portrait) and iPhone (landscape): two header rows** — brand with the data-source
  controls, then the tabs in one row that scrolls sideways and keeps the active tab in view —
  instead of the old stack of a 2-column tab grid and three full-width rows.
- **iPhone (portrait): three short rows** (brand + Live Data + Demo, the run-folder field, tabs),
  and on phones and short windows (≤ 760 wide or ≤ 540 tall) the shell scrolls as a whole: the
  header scrolls away and pins by its tab row (`--header-pin`, measured in `App.jsx`), switching
  tabs opens the new page at its top, and the AI Assistant sizes its chat to the screen below the
  pinned tabs. The header took 267 of the iPhone 17 Pro's 874 px before; it takes 175 at the top
  of a page and 63 once scrolled. In landscape it took 264 of 402 px, leaving 138 for the page;
  now 129, then 63. Card headers wrap their actions under the title on phones (the Displacement
  Directions and Displacement statistics titles no longer run under their controls), `?` help
  opens as a full-width sheet under its trigger instead of running off the right edge, and
  nothing lets a page pan sideways (hidden popovers past the right edge had let iPhone pages
  scroll up to 100px sideways).
- **No cut-off or overlapping text.** An automated audit — every page, with no run and with the
  Demo run loaded (Bond Geometry computed, Dashboard lists expanded), at 22 content sizes from a
  320px phone to 4K — checks for text clipped by its box or ellipsized, text over text, controls
  over controls, cards over cards or the footer, and anything past the screen edge; it ends at
  zero findings. What it fixed: file names under card titles (the χ² log list, the `.rmc6f` name
  in a Model information title cell) and in the loaded-file chips wrap instead of ending in "…"
  (the log list between names); the tolerance ladder keeps every space-group label whole (a
  brick too narrow for its label takes room from the others, with its exact share where all fit,
  and on a phone the bar wraps only if the labels truly cannot share one row); on phones Live Data
  moves up beside the brand so the run-folder field shows "No folder selected" in full; on a
  320px phone the brand keeps its mark only, footer links no longer break inside themselves,
  slider rows shrink the slider, and the PCA camera-axis overlay stacks above the projection
  toggle; and the Atomic Density "Computing KDE…" message is centred instead of printed over the
  previous slab readout in the corner.
- **Stat cards (Model information, Detected SG, Triplet result) never wrap raggedly.** On a phone
  their stats are a two-column grid — aligned, left-aligned, no column dividers — with the cell
  lengths and the tolerance ladder on full rows and the title bar reading like a card header
  (file name at its right end). Before, they wrapped into rows of two or three at random widths,
  each row starting indented behind a stray divider, with the move counters right-aligned under
  left-aligned stats. On wider screens a rail keeps its single row wherever it fits; where it
  does not (iPhone landscape, iPads, 13" MacBooks), `StatRail` measures that and Model
  information breaks into bands — the cell, the atom counts, the move counters — packed into as
  few rows as fit (three on an iPhone in landscape, two on an iPad or a MacBook Air 13"), each row
  left-aligned with dividers only between neighbours. The plan is `ui/statRows.js`
  (unit-tested); bands are marked with `<Stat band>`.
- **Touch and iOS details.** `viewport-fit=cover` with safe-area insets on the header and page
  edges (notch / Dynamic Island in landscape, home indicator); `100dvh` (`--viewport-h`) for the
  fixed shell and every viewport-height budget, so iOS / iPadOS Safari's toolbars no longer hide
  the bottom of a page; taller controls and range thumbs on touch screens (`pointer: coarse`);
  16px text in fields on phones and portrait tablets, so Safari no longer zooms into a field on
  focus; no tap flash or double-tap-zoom delay on controls; an `apple-touch-icon` for the home
  screen.

### Upgrading from 0.5.0

Each 0.5.0 → 0.6.0 example below comes from running both releases on the same input. The inputs are
the committed demo run (`web_app/frontend/public/demo`, GaTa₄Se₈ 250 K, 52 000 atoms), the
maintainer's Mn₃Sn and FeCoSn 199 K total-scattering runs (`data/stog_tests`, not in the
repository), or a synthetic cell where stated. The reasons are in
[Correctness fixes](#correctness-fixes), engine by engine.

#### Numbers that change

| Quantity | Where | 0.5.0 → 0.6.0 | Cause |
| --- | --- | --- | --- |
| Rwp | Dashboard chip; `metrics.rwp` of `/api/plot/*`; `make_plot` | Demo F(Q) 4.569 % → 4.565 %, x-ray G(r) 1.189 % → 1.186 %. A fit with calc = 0.7 × expt: 42.9 % → 30.0 %. The change grows with the amplitude mismatch between the curves. | Divided by the experimental curve, not the calculated one ([parsers](#parsers-and-run-dashboard)) |
| Structure KDE peak height and contour levels | Atomic Density map, both runtimes; `vmax` and `contours` of `/api/kde/slice`; `oriented_kde_slice` | Demo at the default dz 0.08 and bw 0.03 (c, a and Se-only slices): 0 to +23 %, most between +3 and +16 % (c-slice z = 0.55: 115.8 → 134.2). Thick or sparse slabs rise up to +64 % (z = 0.5, dz = 0.2: 49.6 → 80.8; Se at z = 0.3, dz = 0.1, 11 atoms: 220.6 → 362.4). The integrated density is unchanged, and a slab whose atoms all sit on one layer does not move (z = 0.25: 596.4). | The kernel is fitted to the slab's own atoms, without periodic images or the subsample, so it is narrower ([Structure KDE](#structure-kde)) |
| Auto StoG a, b: first shell and fit window | `rmc-autoscale`, `/api/scaling/*`, Auto StoG page; density-limit fits without r₀ | Mn₃Sn (composition given, ρ₀ 0.063049 Å⁻³) at Qmin 1.0: 59438 Q 1.0–27 (a, b) = (0.335, 0.666) → (1.105, −0.102); 300 K 0.331 → 0.752; 500 K 0.272 → 0.936. At Qmin 0.82, 59438 Q 0.82–28: 0.362 → 1.289. | 0.5.0 took the second shell (onset 3.43–3.53 Å) for the first, so its window crossed the first shell ([Auto StoG](#auto-stog)) |
| Auto StoG a, b: Huber weighting | same | Where 0.5.0 already found the first shell, the scale drops: 300 K Q 0.82–28 1.309 → 1.208 (−7.7 %), 500 K 1.567 → 1.397 (−10.8 %), 59438 Q 1.0–29 with `--r0 2.67` 1.199 → 1.011 (−15.6 %). | Rows scaled by √w (the Huber estimator), not by w |
| Auto StoG FZ amplitude and ρ₀ estimate | `--amplitude fz`, `--estimate-rho0` | FeCoSn 199 K (its stog.inp, ⟨b²⟩ = 1.10426): a_fz 1.25116 → 1.1784 (−5.8 %), ρ₀ 0.060692 → 0.057045 Å⁻³ (expert value 0.057329); FeCoSn 100 K ρ₀ 0.0640 → 0.0601 Å⁻³ (+11.7 % → +4.8 %). The density-limit scale barely moves: 1.18293 → 1.18536. | Huber √w |
| Auto StoG with a formula whose ⟨b⟩² disagrees with the stog.inp | CLI `--formula`, API `formula` | FeCoSn stog.inp (⟨b⟩² = 1) with `--formula FeCoSn`: a 1.18874 → 1.18536. The formula's ⟨b²⟩ is no longer used, so `a_fz` and the concordance line (0.5.0: a_fz/a = 0.427, DISCORDANT) are gone unless you pass `--b-sq-avg`. | ⟨b⟩² and ⟨b²⟩ come from one source |
| Automatic low-r enforcement cutoff | the `_rmc` RMCProfile files; the CLI `enforcement:` line | 59438 Q 1.0–27: 3.48 → 2.43 Å; 300 K Q 1.0–27: 3.43 → 2.45 Å; 59438 Q 1.0–29 with `--r0 2.67`: 3.56 Å (above the given r₀) → 2.42 Å. | Cut at the foot of the first shell, never at the onset and never above r₀ |
| PCA `nonGaussianity` | PCA Ellipsoid table; `/api/pca/sites`, `/api/pca/kde` | Demo site 1 (Ta) 4.465 → 4.688, site 2 (Se) 0.456 → 0.398. The per-axis excess kurtosis is unchanged. A symmetric split site now reads negative. | Mean per-axis κ replaced by Mardia's (b₂ − 15)/5 ([PCA ellipsoid](#pca-ellipsoid)) |
| PCA spread of a site at x = ½ of a one-cell-thick box | PCA Ellipsoid, Displacement Directions | Synthetic 8 × 8 × 1 box (c = 10 Å), site at z = ½, σ = 0.08 Å: PC1 rms 4.98 → 0.080 Å. The demo (10 × 10 × 10 cells) is unaffected. | Offsets unwrapped about the site's circular mean |
| Displacement Directions Auto resolution | Resolution "Auto"; `/api/pca/orientation` without `frequency` | Demo sites (1000 copies each): ν = 3 (92 cells) → 2 (42 cells). 300 points: ν 2 → 1; 12 000 points: ν 10 → 9. | `recommended_frequency` floors ([Displacement Directions](#displacement-directions)) |
| Displacement Directions significance | significance strip; payload | Demo site 1 at ν = 10: the strip's "map significance 1.0σ" (the legacy RMS z `significance`) and the peak's "z = 5.4" give way to `peakSignificance` 0.82 σ and `mapSignificance` 0.14 σ. | Calibrated tests replace local z-scores |
| Antipodal-asymmetry null | same | Demo site 1 at ν = 10: `antipodalAsymmetryNull` 0.565 → 0.527 ± 0.018. The asymmetry itself (0.562) is unchanged; z = 1.96, not flagged. | The exact inversion-symmetric null |
| Bond-angle counts, A = C with distinct overlapping windows | Bond Geometry; `/api/triplets`; `rmc-triplets` | Demo Se–Ta–Se, windows 2.3–2.6 / 2.5–2.9 Å: 173 957 → 159 216 angles, mean 109.02° → 109.37°. Equal windows are unchanged (Se–Ta–Se 2.3–2.9 Å: 221 391). | Each physical triplet counts once ([bond angles](#bond-angles)) |
| Bond counts | same | New `uniqueBonds`: demo Ta–Ta–Ta 2.7–3.3 Å keeps `count` 43 624 (B-centred bond vectors, as before) and reports `uniqueBonds` 21 812. | Each physical bond counted once |
| Detected SG ladder | Dashboard Detected SG card | Demo: P1 / Pn / P2 (No. 3, given to a 4-operation set) / P3m / F-43m → P1 / Cm (No. 8) / Cmm2 (No. 35) / P-42₁m (No. 113) / F-43m (No. 216). The full group appears from 0.031 Å (was 0.038 Å). | Closed groups only, standard settings, refined translations ([symmetry finder](#symmetry-finder)) |
| `dispA` | AI Assistant context (`mean_disp_A`, `max_disp_A`) | Isotropic cloud, 0.1 Å per axis: hexagonal cell 0.190 → 0.172 Å, fcc in its 60° rhombohedral cell 0.211 → 0.172 Å. Cubic cells, the demo included, are unchanged. | Mapped through the full cell metric |

#### Files that change

Written files:

- **Classic Auto StoG outputs follow Fortran stog.** `scale.gr` / `<stem>.gr` hold the
  unfiltered g(r), which tends to 1 at large r. `scale_ft.gr` / `<stem>_ft.gr` hold the filtered
  g(r) plus a third column r·[g(r) − 1]. In 0.5.0 they held g − 1 and 4πρ₀r(g − 1). The `_rmc`
  files (Keen G_K, D, F_K) keep their conventions; their values move with a and b (above).
- **The provenance JSON** (`<stem>_provenance.json`) gains `diagnostics.a_fz_rel_se`,
  `a_fz_reliable`, `r_alias_limit` and `rmax_beyond_alias_limit`, `enforcement.source`,
  `provenance.fz_limit` (the Q→0 head fit), `provenance.fit_failure` and
  `provenance.r_alias_limit`, and `rho0_estimate.reason`, `stopped` and `q_first`. When the
  formula's ⟨b⟩² disagrees with the stog.inp value, `diagnostics.a_fz`, `amplitude_concordance`,
  `amplitudes_concordant`, `fk_qmin` and `fk_q0_theory` are absent unless ⟨b²⟩ is given.
- **Output families are written whole or not at all.** `rmc-autoscale`, `/api/scaling/run`,
  `rmc-triplets` and `/api/convert/frac` check every destination before computing and write
  through temporary files renamed into place. A failed run leaves no partial family and never
  replaces an existing file.

Files read differently:

- **stog.inp line 22** follows the classic r₀ rule: `peak_rmin` when the first-peak window starts
  inside the cutoff, else `peak_cutoff`. `2.48 2.65 3.1` now gives r₀ = 2.48 Å (0.5.0: 2.65 Å). A
  value that would leave a fit window narrower than 0.1 Å is ignored, and r₀ is detected.
- **`.rmc6f` atom lines** must be `id element [label]` followed by exactly 7 fields
  (`x y z ref cx cy cz`) or 3 (coordinates only). 0.5.0 misread other layouts silently in both
  runtimes. With one extra trailing field on every demo line, the Flask c-slice at z = 0.25 was
  empty, `/api/pca/sites` found 10 sites instead of 52, and the browser model had 10 basis sites.
  Such lines are now skipped and counted, and a file with none parseable is an error quoting the
  first line. Element names are capitalised in the browser as in Python (SE → Se).
- **A CSV with a stray non-numeric cell** fails with a line-numbered error in the browser. It used
  to plot a silent NaN.

Inputs that now stop before any output is written, with a message naming the problem (HTTP 400, a
thrown worker error or a CLI exit code). 0.5.0 returned a result for most of them, or failed
partway:

- **Auto StoG.** A first shell within ~0.55 Å of `r_cutoff` (Si–O, P–O, B–O and C–O at the
  default 1.0 Å). An unpinned density-limit fit that cannot confirm a first shell: 9 of 56 Mn₃Sn
  (run, Qmin, Qmax) configurations, among them 59438 Q 1.0–28 (0.5.0: a = 0.336) and 55537
  Q 0.82–28 (0.870). A fit with a ≤ 0. Overlapping Q banks, ⟨b²⟩ < ⟨b⟩², a negative or
  non-finite `r_cutoff`, a pinned density-limit window narrower than 0.1 Å. A negative or NaN
  Qmin, a non-finite Qmax, `--scale` NaN or 0, a non-finite `--offset`. An explicit enforcement
  cutoff that is non-finite, negative or at/beyond rmax, and a reversed first-peak window.
- **Auto StoG outputs.** An output that names the input data file or the stog.inp, even with
  `--force` (in `--data` mode, `--out-dir` = the data folder now needs `--out-stem`). Two outputs
  on one file, a directory at an output path, an output folder that cannot be created.
- **`.rmc6f` headers.** A zero, negative or fractional supercell, and a lattice with a
  non-numeric, NaN or missing row, a zero volume or an overflowing volume. `Supercell dimensions:
  0 0 0` used to fold every atom onto the origin.
- **Analysis requests.** A KDE-slice element the file lacks, a zero-spread PCA site's KDE, a
  PCA-KDE volume that captures less than 10⁻⁶ of the density, a negative, NaN or fractional
  orientation smoothing or frequency, a Bond Geometry request above 5×10⁷ angles (app boundaries only; library and CLI
  calls are unrestricted), and a cleared bond-window box.
- **Frac conversion.** A file with no full-layout atom line, and an output that is its own source
  or a directory.
- **`rmc-triplets`.** Destinations that clash, the configuration or a directory as a destination,
  a folder that cannot be created, and a plot format matplotlib cannot write.

#### API changes

| Endpoint | Added | Changed | New statuses |
| --- | --- | --- | --- |
| every `/api/*` route | — | Invalid numbers (NaN, ±∞, text, out of range) are refused by `_number()`; a NaN in a data series is `null`; a non-finite computed result is refused. | 400 where 0.5.0 answered 200 (e.g. `bw=nan` returned `"bw": NaN`, invalid JSON) or 500. 409 when the source file keeps changing during the read. An unknown `/api/*` path is a JSON 404 and a wrong method a JSON 405 with `Allow` (0.5.0: HTML pages). |
| `/api/plot/data`, `/api/plot/metadata` | `chiColumn` (`.log`) | `*_FQn.csv`: title "S(Q) (x-ray)" → "F(Q)", `yLabel` S(Q) → F(Q). `*_PDFpartials.csv`: "PDFpartials", G(r) → "Partial g(r)", g(r). `.log`: "R-value", log(χ), series "R" → "χ² history: X_ray_(R)1", ln(χ²), series "X_ray_(R)1". Classic `scale.gr` / `scale_ft.gr` are labelled g(r). `metrics.rwp` values (above). | An unsupported or unreadable file is a 400 on every plot route (0.5.0 answered some with 500). |
| `/api/structure` | `parseReport`, `parseWarning` | — | 400 for an invalid header or zero parseable atoms. |
| `/api/kde/slice` | `kernel`, `message`, `messageCode`, `warnings`; query `ux`…`vz` (the page's in-plane frame for a custom plane) | `element` is case-insensitive (0.5.0: `element=se` drew an empty map). `bw` is `null` for a declined bandwidth. For a custom normal, `z`/`dz` echo the slider fractions: (1 1 1) at z = 0.5, dz = 0.05 echoed 0.866 / 0.0866 in 0.5.0 and echoes 0.5 / 0.05 now; `depth`/`depthThickness` keep depth units. `z` is clamped to [0, 1] and echoed clamped. A custom (1 1 0)-type map in Flask mode is drawn in the page's frame, rotated 90° from 0.5.0. | 400 for an unknown element (0.5.0: 200 with an all-zero map), an overflowing bandwidth and a non-orthogonal frame. |
| `/api/pca/sites` | per site `axisResolved`, `zeroSpread`, `elementCounts`, `mixed`; top-level `parseWarning` | `nonGaussianity` is Mardia's (above). A zero-spread site has `null` axes, anisotropy, κ and `nonGaussianity`. | 400 for an invalid `.rmc6f`. |
| `/api/pca/kde` | `axisResolved`, `boxHalfWidths`, `elementCounts`, `mixed` | `cubicBox` only sizes `boxHalfWidths`; `halfWidths` is always the per-axis box. | 400 for a zero-spread site or a volume that captures less than 10⁻⁶ of the density (0.5.0: a 200 all-zero volume). |
| `/api/pca/orientation` | `peakSignificance`, `peakPValue`, `peakLocalPValue`, `peakCount`, `peakExpected`, `peakTieCount`, `mapSignificance`, `mapPValue`, `mapChiSquare`, `mapDegreesOfFreedom`, `mapNullSd`, `mapNullSkewness`, `mapExpectedPairs`, `antipodalAsymmetryNullSd`, `antipodalAsymmetryZ`, `antipodalAsymmetrySignificant`, `orientationEffectivePoints`, `orientationAnisotropyNull`, `orientationBinghamStatistic`, `orientationBinghamPValue`, `orientationAnisotropySignificance`, `parseWarning` | `mapPValue` and `mapSignificance` are `null` on maps too sparse to test. `recommendedFrequency` and `antipodalAsymmetryNull` values change (above). The legacy `significance` (an RMS z) and `peakZScore` (a local z) stay but are no longer displayed. | 400 for negative, NaN or non-integer `smoothing` and a non-integer `frequency` (0.5.0: `smoothing=-1` was a silent no-op). |
| `/api/triplets` | `lengths12/23.uniqueBonds`, `parseWarning` | A = C counts with distinct windows (above). | 400 above 5×10⁷ angles, counted before any angle is formed. |
| `/api/scaling/preview`, `/api/scaling/run` | `warnings` (the coefficient warnings the CLI prints); in `diagnostics` `a_fz_rel_se`, `a_fz_reliable`, `r_alias_limit`, `rmax_beyond_alias_limit`, `first_shell_below_r0`; in `provenance` `fit_failure`, `fz_limit`, `r_alias_limit` | `enforce` is one tri-state. The concordance fields can be absent (see Files). JSON booleans are strict: `true`/`false`, `0`/`1` or the words 1/true/yes/on and 0/false/no/off. | 400 for every refusal listed above, a malformed `.inp` (with its own parse error), a non-boolean flag and a body nested too deeply. 404 when a stog.inp names a folder as its data file. 409 when an output exists and `force` is false (unchanged). |
| `/api/convert/frac` | `parseWarning` | `overwrite: "false"` is false (it counted as true). | 400 for the refusals listed above. |

#### Defaults and safety

- **`python web_app/backend/app.py` listens on `127.0.0.1` with debug off.** 0.5.0 listened on
  every interface (`0.0.0.0`) with Flask debug mode on, which exposed the Werkzeug interactive
  debugger (arbitrary code execution) and the data API to the network. Set
  `RMC_TOOLKITS_HOST=0.0.0.0` to listen on the network, or use Gunicorn or Docker. Set
  `RMC_TOOLKITS_DEBUG=1` for local debugging only; the startup line warns when debug meets a
  network address. `PORT` still wins over `RMC_TOOLKITS_PORT`, and a malformed value stops the
  server with an error naming the variable (`server_settings()`).
- **The Docker image no longer bundles `data/`.** The `Dockerfile` used to `COPY data ./data`,
  which failed on a clean clone (the folder is gitignored) and baked private runs into the image
  where it existed. `.dockerignore` now excludes `data/`, `node_modules`, virtualenvs and caches,
  and the image has an empty `/app/data` to mount runs on
  (`docker run -v /path/to/runs:/app/data …`).
- **PCA isosurface level: 25 % → 50 %**, the ellipsoid's probability.
- **Automatic low-r enforcement** cuts at the first-shell foot (numbers above). It was already
  automatic in 0.5.0 whenever no cutoff was given (CLI `--data`, the API without a stog.inp, a
  blank Cutoff on the page), although 0.5.0's `--help` said it was off in `--data` mode. A
  stog.inp or `--enforce-cutoff` still wins, and `--no-enforce` turns it off.
- **Displacement Directions** keeps ν = 10 as the default; only Auto moved (above).
- **`--estimate-rho0`** seeds ρ₀ = 0.05 Å⁻³ when the input gives no density.
- **Live Data in Flask mode** reloads the analysis pages in place when the `.rmc6f` changes. With
  Live Data off, press Load to pick up a newly saved configuration.
- **The static-mode workers refuse what the Flask routes refuse** (see
  [Robustness and parity](#robustness-and-parity)), so a malformed request errors in both runtimes.

### Highlights

- **Physics fixes with visible consequences.** The dashboard Rwp was divided by the calculated
  curve, reading 43 % high for calc = 0.7 × expt. PCA sites at x = ½ in a one-cell-thick box were
  torn in two, giving an rms of half the box edge (5.0 Å in a 10 Å box) instead of 0.08 Å. Auto
  StoG fitted its density-limit window across the first shell of most oxides (synthetic SrTiO₃:
  a = 5.28 for a true 10; ReO₃: a = −3.2), and its automatic low-r cleanup cut 6–9 % of the first
  shell away. The Detected SG ladder named operation sets that are not groups ("P2 No. 3" for a
  4-operation set on the bundled demo).
- **One answer in both runtimes.** The browser ports are held to committed Python goldens: Auto
  StoG to round-off (≤ 5×10⁻¹³ on every real run), the Structure KDE to 10⁻⁶ of the peak
  (measured ≤ 2×10⁻¹²) on slabs within its 6000-row fit cap, bond angles to exact integer
  histograms, and displacement directions to 10⁻⁹ relative. The cap counts the periodic-image
  rows too, so a few thousand atoms reach it at larger bandwidths. Above it the two runtimes sum
  different unbiased subsamples, and their maps differ by several percent of the peak (5.6–8.5 %
  on over-cap slabs of the bundled 52 000-atom demo run); the Flask/SciPy map is the reference.
  The shared `.rmc6f` parser and the plot readers have goldens too. PCA is the exception: each
  engine is pinned to its own reference, with no shared golden ([ALGORITHMS.md](ALGORITHMS.md)).
- **Auto StoG files match classic stog.** The classic-named files hold the Fortran stog
  functions, g(r) and r·[g(r) − 1], and the stog.inp r₀ rule is the classic one. Automatic
  low-r enforcement cuts at the foot of the first shell, and a fit with a non-physical scale is
  never written.
- **Statistics you can quote.** Displacement Directions reports calibrated significances instead
  of local z-scores. PCA non-Gaussianity is rotation-invariant (Mardia). Bond Geometry reports
  physical bond counts and counts each same-element triplet once.
- **Inputs are validated, not guessed.** Both runtimes share one validated `.rmc6f` atom-line
  grammar, with a parse report against the header's atom count; this fixes the zero-atom
  static-mode issue. Every Flask numeric parameter is range-checked (bad input is a 400, never a
  200 with invalid JSON), and parsed-file caches never serve a torn read.
- **Stable-release tooling.** Dependency floors are declared, and CI runs a Python
  3.9/3.11/3.13 matrix. Tests on the committed demo run now exercise real RMCProfile files in CI,
  where the real-data tests skip.

### Correctness fixes

#### Auto StoG

*`scaling.py`, `scaling_cli.py`, `autoScale.js`, `rmc-autoscale`*

- **First-shell detection finds the first shell.** `detect_first_peak_onset` took the flank of
  the strongest |g| feature, which skipped weak or inverted first shells: SrTiO₃'s Ti–O read at
  2.6 Å instead of 1.85 Å, and Mn₃Sn 59438 (Q 1.0–28) at 3.49 Å instead of ~2.74 Å. It now returns
  the smallest-r shell of either sign that stands out of the ripple field below it. A flank that
  reaches the search start is not a shell.
- **The density-limit window is placed from the data.** Without r₀ or `r_fit_max`, the blind
  [1.2, 2.2] Å window used to sit across the first shell. Two trial fits now propose shell onsets.
  The smallest onset is refitted on [lo, onset − 0.25] Å and accepted only when that refit
  re-detects it within 0.15 Å, so a ripple lobe is dropped (at most 4 refits). A refit with
  a ≤ 0 raises. Every fit returned on the four Mn₃Sn runs (Qmin 0.82/1.0 × Qmax 24–30) has a > 0
  and a window top of 2.40–2.50 Å.
- **Automatic low-r enforcement at the first-shell foot.** Enforcement used to zero g(r) up to
  the first-shell onset, deleting 6–9 % of the first-shell pair density. It now cuts at
  min(foot, anchor − 0.25 Å), where the anchor is min(detected onset, given r₀) and the foot is
  `first_shell_foot`. The cut never lands above a given r₀, and the first-shell coordination
  number is preserved to ≤ 0.3 %. On Mn₃Sn 59438 at Q 1.0–27 the cut is at 2.43 Å (expert
  `rmccut` 2.48 Å). A given r₀ that sits above the detected shell is flagged
  `first_shell_below_r0`.
- **Huber IRLS is the Huber estimator.** Both engines multiplied rows by the weight w, which
  minimises Σw²e² (a redescending estimator). Rows are now scaled by √w. On FeCoSn 199 K the
  density and Faber-Ziman amplitudes agree to 0.6 % (5.8 % before), and the ρ₀ self-consistency
  estimate lands 0.5 % from the expert density (5.9 % before). The √w-only table and the
  0.5.0 → 0.6.0 table are in [auto-stog.md](algorithms/auto-stog.md) (Step 6).
- **The page fits what the CLI fits.** The JS loop filtered with S(0) = 0 instead of the
  composition target, so page and CLI (a, b) differed by up to 2.0 % on Mn₃Sn. They now agree to
  round-off, and so does the iterated ρ₀.
- **ρ₀ self-consistency.** `estimate_rho0` confines ρ₀ to 0.005–0.25 Å⁻³ and accepts a root only
  where the density limit holds. With r₀ pinned at 2.6 Å, Mn₃Sn 300 K had "converged" at
  0.428 Å⁻³, about 50 g/cm³. Every
  non-converged exit carries a `reason`, and the CLI and page name the trial density that could
  not be fitted. `extrapolated` is judged on the first measured Q.
- **⟨b⟩² and ⟨b²⟩ come from one source.** A composition's Sears ⟨b²⟩ was paired with an
  unrelated ⟨b⟩² (for example 1 for normalised x-ray data). That invented an S(0) target (+0.55
  for FeCoSn) and a 60 %-low "converged" density. S(0) > 0 is now rejected.
- **Q order and bank overlap.** Descending-Q files are sorted. Before, they gave a negated G(r),
  and a reversed synthetic fitted a = −10.1 for a true +9.97. Duplicate or overlapping Q
  (concatenated banks) raise.
- **Despiking and σ.** `--despike` runs once, so the fit and the written files share one point
  set (2136 vs 2308 points on 59438 before) and `n_despiked` is the true count. The page applies
  the CLI's σ-column guard, and an unconverged JS run reports `max_iter` iterations.
- **Faber-Ziman conditioning.** New fields `a_fz_rel_se` (Huber sandwich error) and
  `a_fz_reliable`. A reliable a_fz is described as necessary, not sufficient, and 59438's a_fz is
  flagged unreliable.
- **Aliasing limit.** `r_alias_limit` is π/max ΔQ of the transformed grid, and an r_max beyond it
  is flagged. Despike gaps lower it to 24 Å on 59438.
- **No sliver fit windows.** A stog.inp line-22 cutoff pins r₀ only when the default window it
  leaves is at least 0.1 Å wide (`MIN_AUTO_WINDOW`, the automatic placement's floor); otherwise
  r₀ is detected. A cutoff just above r_cutoff + 0.45 Å used to pin a 1–6 point window: on
  FeCoSn 199 K, `1.46 0 0` gave a = 1.001 and `1.0 0 0` with r_cutoff 0.5 gave a = 0.672 (1.185
  with the shipped file), both "converged" with the density limit "satisfied". A pinned
  density-limit window (r₀ / r_fit_max) narrower than 0.1 Å now raises. `r_cutoff` must be finite
  and ≥ 0 in both engines (a negative value returned a "converged" a = 0.026 and was written).
- **The input is never an output.** In `--data` mode the default stem is the data file's, so
  `--out-dir` pointing at the data folder made `<stem>.sq` the measured file itself, and
  `--force` (which the refusal suggested) replaced it with the scaled, cropped S(Q); a rerun then
  silently fitted already-scaled data (a = 1).
- **Outputs are written whole or not at all.** A stog.inp declaring the FK(Q) name as `ft.dat`
  used to exit 0 with the RMCProfile input replaced by the filter correction, and a declared
  `sub/rmc.gr` or a directory at `ft.dat` failed after five or seven files were written. Every
  destination is now checked before any computation (names compared case-insensitively), the
  family goes through temporary files, and declared subfolders are created. An API `outDir` that
  is a file is a 400, not a 500 after the fit.
- **Input checks in every entry point.** `--scale nan` exited 0 with nine all-NaN RMCProfile
  files. A NaN Qmax passed Python's `qmax <= qmin`. `--enforce-cutoff nan` was reported as
  applied, and a cutoff of 1000 with rmax 50 replaced the whole G(r). `ScalingConfig` /
  `makeConfig` and `validate_enforcement` / `validateEnforcement` (CLI, API, page worker) now
  refuse these. A given r₀ with no detected shell no longer crashes with a TypeError.
- **stog readers.** Both engines read CR-only line endings, a BOM, bad bytes, Fortran D exponents
  and ASCII-only digits the same way.

#### Transforms

*`transforms.py` and its JS port*

- The Lorch low-Q correction basis is cancellation-free. Grid points 10⁻⁹–10⁻⁶ Å from π/Qmax
  got coefficients off by up to ~80× or with flipped sign. The unwindowed moments use a Taylor
  series at small Q₀r, where they were off by 100 %.
- The Fourier filter accepts an r grid that starts at 0. Before, every output Q was NaN;
  `g_filtered(0)` is now the continuous extension (`gpdf_slope_at_zero`).
- The public transforms raise on a grid that is not strictly increasing.

#### Structure KDE

*`kde.py`, `localKdeWorker.js`, `gpuKde.js`*

- **The browser draws the SciPy kernel.** The worker's 10⁻⁸ ridge and its isotropic fallback are
  gone. Both runtimes evaluate H = bw²·C exactly and decline the same slabs with the same
  `message`. Before, on the GaNb₄Se₈ Ga layer at bw ≤ 0.015, the browser map differed from SciPy's
  by 36–55 % of the peak.
- **Kernel from the source atoms.** C is the covariance of one row per slab atom. Periodic images
  and the 6000-row subsample enter only the density sum, so the thickness and bandwidth sliders
  no longer reshape the kernel. The kernel is narrower, so peaks rise in both runtimes (sizes in
  [Upgrading](#numbers-that-change)).
- **One slab test**, `|d − z_c| ≤ dz/2 + 10⁻⁹`, in `kde.py`, the worker and the Slab-In-Cell
  highlight.
- **Contours.** Flask's log mode no longer drops contours when the peak is below 1. A map whose
  kernel misses every grid node is flagged `unresolved` and is neither contoured nor painted.
- **SciPy.** `_FixedCovarianceKDE` drives SciPy 1.8 through 1.18; `/api/kde/slice` used to answer
  500 on SciPy < 1.10.
- **Honest labels.** The custom slice is labelled as the Miller plane (h k l) it is. The map
  prints the slab thickness and the kernel σ in Å and warns on sub-grid kernels.
- **Custom (hkl) slices are drawn in one frame.** For a custom plane the page now sends its
  in-plane frame to `/api/kde/slice` (`ux`…`vz`), and the route draws the map in it
  (`_custom_slice_frame`, validated orthogonal to the normal). Flask used to pick its own axes, so
  the default (1 1 0) map — and (1 0 1) — came back rotated 90° against the browser map and against
  the page's own Slab In Cell panel, letterboxed in a panel sized for the other orientation. The
  densities were always the same; the parity golden now includes (1 1 0) and (1 0 1) demo slices.

#### PCA ellipsoid

*`pca_kde.py`, `pcaKde.js`, `pcaCrystalFrame.js`*

- Offsets are unwrapped about each site's circular mean, not folded about zero. A site at x = ½
  in a one-cell-thick box, or near x = 1 in a two-cell box, was torn in two, making U 10³–10⁴×
  too large. Displacement Directions inherited the fix.
- The browser ellipsoid scale uses the exact χ²₃ quantile. The old table was up to +3.5 % off and
  had a wrong 0.6827 node.
- The volume, mass, iso levels and PC walls are always sampled on the per-axis box. On planar
  clouds the captured mass had read 0 % or ~2×10⁵ %. The shell no longer paints clamped box-face
  density, and the crystal-frame walls are line-integral marginals (no moiré; error 6–22 % → 2–3 %
  at grid 40).
- Mixed-occupancy sites are labelled by their majority species. Element pooling selects atoms by
  their own element.
- Both engines refuse a volume that captures less than 10⁻⁶ of the density (a bandwidth far below
  the node spacing, or an extent of 10⁶); it used to be an all-zero volume with no warning.
- The Cartesian → CIF U_ij formula in the docs is corrected, as is the claim that the
  20 000-point cap binds only pooled clouds (it also binds single sites in boxes of ≥ 28 cells per
  edge).

#### Displacement Directions

*`orientation.py`, `orientation.js`*

- **Calibrated significances.** The peak readout is `peakSignificance`: the exact Poisson tail
  of the peak cell, Šidák-corrected over all cells. The map test compares Pearson's X² with a
  gamma matched to its exact isotropic moments, in `mapSignificance`. It replaces an RMS z
  printed as "σ": a hemisphere-only cloud used to read "1.4 σ" and now reads > 10 σ. The map test
  is withheld below 0.1 expected coincident pairs.
- **Antipodal asymmetry.** The null is now the exact inversion-symmetric one, conditional on
  each pair's total, and carries an engine flag at 3 SD. The old √(C/πN) floor could exceed 1,
  so the red flag could never fire at ν = 10.
- **Anisotropy.** The isotropic expectation is 9/√(10πN_eff), and a Bingham test is added.
- **Ties.** Exact Voronoi ties are broken centrosymmetrically (an exactly centrosymmetric cloud
  read an antipodal asymmetry of 1.000 at odd ν),
  tied peaks resolve by one 10⁻⁹ rule, and neighbour cycles start at their smallest index.
  Before, 185 of 1002 neighbour rows differed between engines at ν = 10.
- `recommended_frequency` returns the largest ν with ≥ 12 points per cell, as its docstring
  promised.
- Non-finite rows and options are rejected by name. Before, Python failed with a LAPACK error
  and JS returned NaN axes.

#### Bond angles

*`triplets.py`, `triplets.js`, `rmc-triplets`*

- **Bounded work.** Angles stream into the histogram, and both app boundaries refuse a request
  whose exact angle count exceeds 5×10⁷ before forming any. The 15 Å cap had allowed ~10⁹-angle
  requests: tens of GB in Flask, `RangeError` in the worker.
- **One count per physical triplet.** For A = C with distinct windows the old ordered rule
  counted overlap triplets twice.
- **Ideal geometries are deterministic.** An angle within 10⁻⁹° of a bin edge bins on it, and a
  bond within 10⁻⁹ Å of a window bound is inside. A shell typed at its exact distance used to
  lose up to 60 % of its bonds.
- **Physical bond counts.** `uniqueBonds` counts each bond once; the directed count is 2× when
  A = B.
- A cleared window box is an error. It used to become rmin = 0 silently: Nb–Nb–Nb 3.5–4.6 Å
  became 0–4.6 Å.
- `rmc-triplets <run folder>` picks the same `.rmc6f` as the app.
- The docs correct the sin θ rationale and give the factor to RMCProfile's `norm/sin(theta)`.
  Angle totals match RMCProfile's own TRIPLETS output for five triplets on the 5 K run.
- **`rmc-triplets` checks its destinations before computing.** `--dump-angles` equal to
  `--output` replaced the histogram with exit 0, and an unsupported `--plot` format printed a
  traceback after the CSV was written. `--version` added.

#### Parsers and run dashboard

*`parsers.py`, `plots.py`, `browserData.js`, `rmc6f.js`*

- **Rwp** is normalised by the experimental column. The column roles come from the header when
  it names them, else from RMCProfile's (x, calc, expt) order.
- **`.rmc6f` atom lines.** Both runtimes read the fields from the end of the line, so one extra
  trailing field shifted every coordinate by one column, silently (examples under
  [Files that change](#files-that-change)). Both now share one anchored, validated grammar. It
  reads D exponents, any spelling of the `Atoms` marker, bare-CR files and coords-only lines. It
  skips and counts unparseable and non-finite lines and reports "k of n atom lines unparsed".
- **`.rmc6f` header.** `read_cell_vectors` and `readRmc6fCellVectors` read the header numbers
  like the atom lines (D exponents) and validate the supercell and the lattice, with the same
  error text in both runtimes. A NaN, collinear or 10³⁰⁰ lattice used to give 200s with NaN
  cells, an anisotropy of 10¹⁴ or all-zero angle counts, and a D-exponent lattice was a raw float
  error.
- **Structure files.** 0-byte or marker-less candidates are skipped; an empty `.rmc6f` used to
  hide six valid configurations in `data/250K_try1/supercell`. `read_structure` pairs Frac/rmc6f
  files by stem. Before, it could fold a 5×10×10 Frac file with a 10×10×10 supercell.
- **χ² history.** `.log` reads check the header's column count, drop a half-written last line and
  keep NaN rows as gaps. Static mode charts one run's logs, not every run spliced together. The
  curve is named by its column (`χ² history: X_ray_(R)1`): it is one fit term, not a total.
- **Labels.** `*_FQn.csv` is labelled F(Q), `*_PDFpartials.csv` g(r); PNG and Flask use the same
  strings.
- **dispA** goes through the full cell metric. It was +10 % in hexagonal and +23 % in
  rhombohedral cells.

#### Symmetry finder

*`symmetry.js`, `spaceGroupSymbol.js`, `wyckoff.js`*

- **Closed groups only.** A space-group number is given only to an operation set that is closed
  under composition. A set that fails is labelled "not a group" with no number (1 of 2065 rungs in
  random sweeps). Each operation's translation is least-squares refined, so the answer no longer
  depends on atom order; the demo's full-group residual fell from 0.038 to 0.031 Å.
- **Standard settings.** A symbol is shown only when it is a tabulated symbol of the detected
  class and centring, in a standard setting found from the group's own elements (any axis order,
  centred, primitive or rhombohedral cells, supercells reduced). Otherwise the card shows the
  crystal class or a "≥" lower bound, with no number. A centred group given on a primitive cell
  is named in its centred cell (C2 used to read P2 No. 3, I4 P4 No. 75, R3 P3 No. 143).
- **Lattice rotations** are tested by Cartesian strain on the τ scale. A c/a = 1.004 perovskite
  reads `P4/mmm` below 0.016 Å and `Pm-3m` above it.
- **I2₁2₁2₁ and I2₁3 are told apart** from I222 and I23. The 0.5.0 entry listed that as a known
  limitation.
- **Wyckoff letters** are read in the naming cell and paired with that cell's multiplicity. The
  assistant's labels follow.
- **Correction to the 0.5.0 entry.** The table holds 1731 positions, not 1724. The seven
  P222₁ (#17) and I2₁2₁2₁ (#24) special positions the 0.5.0 check "rejected" are correct. The
  rejection came from test fixtures that described those two groups about a shifted origin. At
  the ITA origin, P222₁'s 2a site `x,0,0` lies on the 2-fold axis along a and expands to 2
  points. All seven are restored.
- **Bounded cost.** Bases over 2000 sites, or boxes over 384 candidate operations, read
  "not analysed" instead of freezing the page. An unanalysable lattice reads "undetermined"
  instead of P1 No. 1.

#### Flask API

*`web_app/backend/app.py`*

- **One validator for every numeric parameter** (`_number()`). NaN, ±∞, text, lists and
  out-of-range values are a 400 that names the parameter. Before, they could produce a 200 with
  invalid JSON or an all-zero map. A computed result that comes out NaN/∞ is also a 400,
  including a KDE bandwidth whose kernel overflows and a scaling fit that overflows.
- **Parsed-file caches** are keyed on the full file signature (`st_mtime_ns`, `st_ctime_ns`,
  `st_size`, `st_ino`), never on a whole-second mtime. A file rewritten over sshfs or `scp -p`
  used to be served stale until the server restarted. A parse of a file that changed while it
  was being read is never cached; a file that keeps changing is a 409.
- **Live Data in Flask mode** re-checks the `.rmc6f` signature on every poll and on Load. The
  analysis pages then reload in place, keeping the picks that still apply, so a page never mixes
  two configurations. Only Bond Geometry's computed distribution is dropped. Every Three.js view
  releases its WebGL context on teardown.
- **Request edges are 4xx, never 500 or a silent 200.** The static-asset rule matched `/api/*`
  before the SPA route, so unknown paths and wrong methods got Flask's HTML pages. A malformed
  `.inp` with the default `kind: 'auto'` was re-read as S(Q) data (preview: "data mode requires
  qmin and qmax"). JSON booleans were read loosely (`"maybe"` read as false, an object as
  true, `inspect: "false"` entered inspect mode). A stog.inp whose data file is a folder and a
  JSON body nested thousands deep were 500s. `/api/kde/slice` matched `element` case-sensitively
  and drew an all-zero map for an element the file lacks. `/api/convert/frac` wrote a
  header-only Frac file with a 200 for a file with no full-layout atom line, and answered 500
  for an output that is its own source or a directory. The new statuses are in
  [API changes](#api-changes).
- `/api/scaling/*` parse `enforce` once, as a tri-state. They return the coefficient warnings
  the CLI prints and never mutate a cached result.
- [REFERENCE.md](REFERENCE.md) documents all 15 routes with their ranges and caps, and a
  contract test keeps it complete.

### Robustness and parity

- Python ↔ browser goldens: Auto StoG (`autoscale_fixture.json`, round-off), Structure KDE
  (`kde_parity_fixture.json`), bond angles (`triplets_fixture.json`), displacement directions
  (values shared between `test_orientation_fixes.py` and `orientationFixes.test.js`), the
  `.rmc6f` parser (`rmc6f_coords_only_fixture.json`) and the plot readers
  (`plot_parity_fixture.json`, including neutron, Bragg and EXAFS layouts). The WebGPU shader is
  checked by a float32 emulation.
- Non-finite input is handled one way everywhere: `.rmc6f` lines with NaN/∞ coordinates are
  skipped, counted and reported in both runtimes (no batched `LinAlgError`). NaN in a data series
  is `null`. A NaN result is a 400.
- Both runtimes resolve a flat run folder to the same configuration (`find_run_configuration`,
  `chooseStructureFile`), with code-point tie-breaks. A picked folder with subfolders can differ:
  the browser also searches the subfolders (and falls back to the first usable `.rmc6f` by full
  path), while the server looks only at the folder itself.
- The package root exports the 0.6.0 engine API (`first_shell_foot`, `auto_enforcement_cutoff`,
  `fz_limit_fit`, `alias_limit`, `bond_angle_summary_from_file`, …).
- The static-mode workers refuse what the Flask routes refuse (`workers/requestGuards.js`): the PCA
  worker's orientation `frequency` must be an integer and `smoothing` an integer in [0, 64] (10⁹
  passes used to pin the worker), a KDE or orientation result holding NaN/∞ is an error with
  `_strict_result_response`'s message instead of a posted NaN volume, and an unknown `kind` is an
  error instead of a silent KDE. The Auto StoG worker requires a finite, non-zero `a` and a finite
  `b` in manual mode (a NaN `a` used to post `ok: true` with all-NaN curves) and refuses a
  non-finite fit. The structure worker reads `maxPoints` as `/api/structure` does. Every worker
  answers a null message with an error, so the caller's promise always settles.
- Coordinates-only site reconstruction on the PCA page is 3.5–6.5× faster. The `.rmc6f`
  classifier's fast path parses 52 000 atoms in ~0.2 s.

### Tooling

- Dependency floors `numpy>=1.22`, `scipy>=1.8`, `matplotlib>=3.6` and `contourpy>=1.0.7`
  (`kde.py` imports contourpy directly). The suite passes on the floors (Python 3.9) and on the
  newest releases (Python 3.13).
- CI runs the Python suite on 3.9 with the floors pinned exactly, and on 3.11 and 3.13 with the
  latest releases. The frontend job runs lint, `npm test` and the build.
- The committed demo run (`web_app/frontend/public/demo/GTS_250K.*`) backs tests that run in CI:
  `tests/test_parsers_demo_run.py`, `__tests__/demoRun.test.js`, the KDE golden and the plot
  parity golden. The GaNb₄Se₈ and `stog_tests` real-data tests still skip without `data/`, and
  the full Mn₃Sn sweep is opt-in (`RMC_TOOLKITS_FULL_SWEEP=1`).
- Version 0.6.0 (`0.6.0.dev0` until it ships; Beta classifier). The Auto StoG tab stays hidden
  in the shipped build (`SHOW_AUTO_STOG = false`); the engine, CLI and API are supported.

### Deferred beyond 0.6.0

The maintainer decisions the audit considered and did not take for 0.6.0 (a physical KDE kernel,
plotting every χ² column, a Python symmetry finder, Auto as the default orientation resolution,
AXIS_RESOLUTION_SIGMAS = 2, automatic FZ fallback for degenerate density limits, bounded A = C
pairing work, and more) are listed per engine in [ROADMAP.md](ROADMAP.md#candidates-after-060).

## v0.5.0 — 2026-08-14

Two new analysis pages, correct space-group naming, and the math reference.

- **Bond Geometry** — bond-angle (triplet) distributions in the app, in both runtimes, on the new
  `rmc_toolkits.triplets` engine and its `rmc-triplets` CLI.
- **Displacement Directions** — the hex-binned orientation sphere.
- **Auto StoG** — composition-first absolute scaling, ρ₀ self-consistency and the `rmc-autoscale`
  CLI. The engine and CLI ship; the tab stays hidden until the page behaves.
- **Symmetry** — non-symmorphic groups named as themselves against all 230, with Wyckoff letters
  per orbit, plus the first tests for the symmetry code.
- **[docs/ALGORITHMS.md](ALGORITHMS.md)** — a code-anchored account of every operation each page
  performs, with seven per-page derivations.
- **PCA Ellipsoid** — principal axes reported in the crystallographic frame.
- Fixes: Rwp reports "unavailable" rather than a fake perfect fit; the ρ₀ fixture tolerance is
  per dataset.

**"Local Geometry" renamed to "Bond Geometry"** (2026-08-14) — the tab shipped in this release
under the name *Local Geometry*, which is generic: in the PDF/total-scattering sense every
analysis page here is local structure, so the name didn't say what distinguishes this one. The
page is bond angles, bond lengths, coordination and detected bonds — bond geometry — and the tab,
component (`BondGeometryPage.jsx`/`.css`), algorithm reference
([algorithms/bond-geometry.md](algorithms/bond-geometry.md)) and every cross-reference now carry
that name. The internal page key stays `geometry`, and the engine keeps its RMCProfile-rooted
`triplets` naming (`rmc_toolkits.triplets`, `rmc-triplets`, `/api/triplets`,
`workers/triplets.js`) — the rename is presentation-layer only. Entries below predating the
rename keep their original titles; the tag was re-cut so v0.5.0 ships the tab under its final
name.

**Local Geometry: folded-cell canvas overflowed the page on retina displays** (2026-08-14) —
`FoldedCellPanel` resized its renderer with `setSize(width, height, false)`, which skips styling
the canvas, and `.pca-structure` gives the canvas no CSS size of its own — so the canvas laid out
at its attribute size, CSS size × devicePixelRatio, twice its panel on a 2× display, and the
oversized canvas wrecked the whole tab's layout. Plain `setSize(width, height)` now lets three.js
style the canvas, as `SiteStructurePanel` always has. Also fixed in passing: the panel's *Reset
view* called `controls.reset()` with no saved state, snapping the camera to its construction
position at the origin (inside the cell); the rebuild now ends with `controls.saveState()` so
reset restores the framed view of the current data. The v0.5.0 tag was re-cut onto this fix's
merge commit; the original 2026-08-13 tag, which pointed at the release-notes commit and so
predated both this fix and the algorithm reference, was deleted before anything referenced it.

**Local Geometry algorithm reference** (2026-08-13) — new
[algorithms/bond-geometry.md](algorithms/bond-geometry.md) (created as `local-geometry.md`,
renamed with the tab), completing the per-page set: every
analysis page now has its own code-anchored derivation document. Engine section: inputs and the
half-up bin-count rule, folding + image bookkeeping, the linked-cell search with the
strictly-more-than-$(q-1)$-thicknesses covering argument, bond admission (inclusive windows,
zero-length exclusion, self-image rule), the two pairing rules, the three normalizations (with
why the bin-integral `sinth` stays finite at 0°/180°), the summary payload, the app-boundary
caps, and the CLI. Page section: seeding, the epoch guard, the `fit` plot variant, the
partial-g(r) helper rules (curves by bond type, guides by window role, palette-matched colors),
and the folded-cell bond view — including the caveat that its sticks are average-structure
bonds, not instantaneous ones. Parity is quoted from the golden-fixture tests (exact integer
counts; 1e-9 / 1e-7 / 1e-5 float tiers) rather than asserted. `ALGORITHMS.md` index row and the
README's bond-angle derivation link now point at it.

**Auto StoG: `test_estimate_rho0_near_hand_value` failed on the 100 K dataset** (2026-08-13) —
the x-ray fixture tests select the first available FeCoSn run (`100K` or `199K`), on the
stated grounds that the two share a stog parameterization. True for every assertion in
that class except one: how closely `estimate_rho0` reproduces the expert's hand density is
a property of the measured S(Q), not of the parameterization. The 4.7 % figure recorded
here and in the plan was measured on **199 K**; on **100 K** the estimate lands at 0.0640 Å⁻³,
11.7 % from the hand 0.057329 — over the hardcoded 10 % bound. Anyone holding only the
100 K run saw a red suite, and CI never caught it because both datasets are gitignored and
the test skips.

- Not an algorithm fault: on 100 K the estimate converges to 0.064033 ± 0.00002 from seeds
  spanning 0.02–0.20, concordance → 1 to 10⁻⁵ in 2–3 iterations, `extrapolated=False`. It
  is the correct root of `a_fz/a_density(ρ₀) = 1` for that data. `scaling.py` is unchanged.
- The tolerance is now per run (199 K 10 %, 100 K 13 %) rather than one bound covering
  both — a shared ~13 % bound would stop the test noticing a regression on 199 K.
- New `test_estimate_rho0_is_seed_independent` asserts the part that does *not* depend on
  which temperature is on disk: seeds spanning a 10× range must agree to <1 %. A
  seed-dependent answer would mean the fixed-point iteration is terminating on its seed
  rather than on the data, which no tolerance on the hand value would catch.
- Both figures are now documented (`docs/algorithms/auto-stog.md` Step 10, plan §ρ₀); the
  100 K entry previously recorded scale accuracy but never ρ₀, since `estimate_rho0`
  postdates it.

**Non-symmorphic space groups are named correctly; Wyckoff letters for all 230** (2026-08-06) —
the *Detected SG* panel reported the **symmorphic parent** of every group with a screw axis or
glide plane. The operations were always right; only the naming step threw the translation parts
away, so diamond came back as Fm-3m rather than Fd-3m, a Pnma structure as Pmmm, an a⁰a⁰c⁻ tilted
perovskite as I4/mmm rather than I4/mcm, and hcp as P6/mmm. Since octahedral-tilt phases are almost
all non-symmorphic (I4/mcm, Pnma, R-3c, Imma), this hit the common RMC case.

- New `spaceGroupSymbol.js`: splits each operation into its **intrinsic translation**,
  `(1/n)·Σ Rᵏ·t`, the projection of `t` onto the invariant subspace of `R`. It is
  origin-independent, which is what makes it a reliable element label — non-zero along the axis
  means a screw `n_m`, non-zero in the mirror plane means a glide `a/b/c/n/d/e`. The symbol is
  assembled positionally over each crystal system's symmetry directions. Handedness comes from
  `det[axis, v, R·v]`, so the enantiomorphic pairs (4₁/4₃, 3₁/3₂, 6₁/6₅, 6₂/6₄) stay distinct.
- New `spaceGroupTable.js`: all 230 groups, replacing a 64-number symmorphic-only table. Pre-2002
  spellings (Cmca → Cmce) resolve to the same number.
- The builder **proposes ranked candidates and the table disposes**, rather than committing to one
  spelling. Necessary because H–M plane priority (m > e > a > b > c > n > d) is not always the
  answer: I-centring puts both b- and c-glides in the same plane of I4/mcm, where ITA writes `c`.
  Non-standard axis settings are searched the same way — a structure handed over in the Pbnm
  setting is reported as Pnma.
- **R centering is detected** (obverse *and* reverse). `CENTERING_SETS` had only F/I/A/B/C, so a
  rhombohedral structure reported `P`; R-3m came out as P-3m1 — wrong centering, and order 12
  claimed for a group with 36 operations. The `R3`/`R-3m` rows of the old number table and the
  `trig: 'PR'` entry in `ALLOWED_CENTERING` had been unreachable.
- A symbol is only built from elements the positions can place. A subgroup part-way up the
  tolerance ladder keeps its cubic parent's axes, putting its 3-fold along ⟨111⟩ where the
  hexagonal-axes positions expect [001]; that used to emit the meaningless `P2m`, and now falls
  back to the crystal class (`P3m`). A symbol merely in a non-standard setting (`Pn` for #7) is
  still shown as-is. Naming also runs only for operation sets that close, so the ladder no longer
  pays for it on partial mid-transition sets.
- **Wyckoff letters for all 230 groups**, up from four (216/225/229/221). New `wyckoff.js` +
  `wyckoffTable.js` (1724 positions, packed one string per group). Multiplicity and site symmetry
  alone are often ambiguous — Pm-3m has both 3c and 3d at multiplicity 3 with site symmetry 4/mmm —
  so ties are broken on the position's coordinate form, tested against every member of the orbit
  because the detected representative is whichever atom came first, not the one International
  Tables prints. That fixes cubic SrTiO₃, whose oxygen orbit previously matched both and got no
  letter; it is now 3c. A letter is still withheld when nothing fits uniquely, since the finder
  never shifts the origin to the standard setting.
- Knock-on fix: diamond used to be identified as #225 and was therefore handed **Fm-3m Wyckoff
  letters for an Fd-3m structure**. It now resolves to #227 and gets that group's own letters.
- **First tests for the symmetry code** — there were none, and the one test mentioning a space
  group hard-codes `F-43m` as an LLM-context fixture without ever running the finder. 736 new
  tests across `__tests__/symmetry.test.js` and `__tests__/wyckoff.test.js`, built on **all 230
  space groups** given by ITA generators in `__tests__/fixtures/spaceGroups.js`. Each group is
  asserted three ways: its generator closure reproduces the general-position multiplicity, its
  symbol is built correctly from exact operations, and the finder recovers it end-to-end from a
  structure made by expanding a generic point. Every tabulated Wyckoff position is re-derived from
  the group's own operations at test time, so the committed data cannot drift from the symmetry it
  claims to describe. Plus an integration test over the bundled 52 000-atom demo model
  (GaTa₄Se₈ → F-43m, Ga 4c / Ta 16e / Se 16e ×2), and unit tests for intrinsic translations, glide
  letters, screw handedness, centering, the ladder, and site orbits.
- All reference data was generated, then **kept only where it verified computationally**: the 230
  symbols come from two independent passes that agreed on every entry; each generator set was
  accepted only when its closure size matched (crystal-class order) × (centring multiplicity); each
  Wyckoff position only when its coordinate expanded to exactly the stated multiplicity. That check
  rejected 7 wrong positions — for P222₁ the proposed 2-fold site `x,0,0` lies on no axis and
  expands to 4 points, not 2.
- Known limitation, pinned by test rather than hidden: **I222/I2₁2₁2₁ and I23/I2₁3 cannot be told
  apart.** I centring gives both members of each pair the same element types along the same
  directions (every pure 2-fold has a 2₁ screw half a cell away), and they differ only in where
  those axes sit relative to one another — which needs the origin-aware matching this finder does
  not do. It reports I222 and I23.
**Local Geometry tab — bond angles in the app, both runtimes** (2026-08-06) — the triplets engine
(merged in #36) is now a workspace page: pick A–B–C with **B central**, bound the bond lengths,
Compute. Payload contract is `rmc_toolkits.triplets.bond_angle_summary` (angle histogram with
counts / per-degree density / exact sin-corrected, bond-length histograms per window, coordination
distribution), served identically by the new Flask `/api/triplets` route (lru-cached on
path+mtime+params) and by the shared PCA worker's new `triplets` request — the same
worker-or-API routing as the PCA pages via `useSiteCloud`.

- `workers/triplets.js` is a line-for-line port of the Python engine (linked-cell periodic search
  with explicit image shifts, shared-ends unordered counting, zero-length pairs never bond),
  parity-tested against Python goldens: exact integer histograms, 1e-9 float agreement
  (`triplets_fixture.json`, regenerated by `tests/generate_triplets_fixture.py`).
  `pcaKde.js` now keeps the per-atom element + supercell-fraction list on its parse result so the
  worker answers triplets requests from the same cached parse.
- `LocalGeometryPage.jsx` follows the PCA/Orientation page pattern (top controls bar, `pca-panel`
  cards, InfoBadges): triplet selects seeded per dataset (ends = most abundant element, central =
  next), inclusive windows with an optional distinct B–C window, result chips (bonds, mean length,
  coordination mode share, mean angle), a sin-corrected|density toggle on the hero plot, and a
  window helper that plots the run's partial g(r) (`PDFpartials.csv`, both modes) with dashed
  guides at the current bounds, cropped to the first-shell region. Three equal-width full-height
  columns: angle distribution, partial g(r), folded unit cell. An earlier fourth column held a
  histogram of the bond lengths inside the window — dropped as redundant, since it was the same
  first-shell peak the partial g(r) already plots, only clipped to the window and without the
  surrounding context that makes the guides worth reading. The engine still returns those
  histograms; the result chips report their counts and mean length.
- **The partial-PDF helper carries both bonds of the triplet.** A second curve is drawn whenever
  A–B and B–C are different bond *types* — Ga–Ta–Se brackets a Ga-Ta and a Ta-Se shell, Se–Ta–Se
  only ever has one, and the pair is looked up in either label order. That is independent of the
  **distinct B–C** switch, which governs the *windows* rather than the curves: with it off one pair
  of guides covers both bonds, with it on each window gets its own pair, labelled `A–B`/`B–C` by
  role — the pair names would be identical for a same-element triplet, whose two windows still sit
  on one shared curve. The first-shell crop widens to clear whichever window reaches further.
  Two windows also take the two curve colors (`plotPalette.js`, new module: a component file that
  also exports constants breaks Fast Refresh). Guides consume no palette slot, so shell N is
  `PLOT_PALETTE[N]` — each window is drawn in the color of the shell it brackets, and a lone window
  keeps the neutral guide stroke since there is nothing to tell it apart from.
- **The 3D panel is the folded unit cell, not site ellipsoids** (`FoldedCellPanel.jsx`): every atom
  folded back into one cell as a point cloud, the same view the Atomic Density page shows, so the
  spread around a site is the measured thermal cloud rather than a fitted ellipsoid. The detected
  bonds are drawn over it as thin transparent lines (`LineSegments`, opacity 0.45) instead of
  opaque cylinders, so a full coordination network reads as a framework without hiding the cloud.
  Bond matching is unchanged — average site pairs inside the window, periodic images included.
  `SiteStructurePanel` keeps its ellipsoid view for the PCA and Orientation pages and gains a
  `title` prop rather than being renamed for one caller.
- **`InteractivePlot` gains an opt-in `fit` variant** that takes its `viewBox` from the rendered box
  (ResizeObserver, rounded to whole pixels) instead of a fixed 8:5 or 1440×320 aspect. The fixed
  aspect letterboxed inside a card of any other shape — the two plots here filled 52–57% of their
  cards, the rest dead space. Under `fit` one user unit is one CSS pixel, so margins keep their
  physical size and tick density follows the box (~1 y tick / 70px, 1 x tick / 95px, clamped);
  measured fill is now 100%. Opt-in because it needs a caller that gives the plot a definite
  height: the fixed variants are untouched and every other page renders exactly as before.
- Both plots are handed **identically sized boxes** so the same data does not read at two aspect
  ratios side by side. Equal columns are not enough on their own — it is the chrome above each
  stage that differs, so the card header is pinned to one height (one title wraps, the other does
  not) and two legend rows are reserved in both (only the partial PDF's legend wraps, once the
  split adds a second pair of guides). Measured: 492×363 in both, viewBox matching.
- The bond-angle normalization note is rewritten. It previously led with the correction rather than
  the thing being corrected; it now defines **density** first, states plainly that random bonds do
  *not* give a flat line there (fewer ways to form an angle near 0°/180° than near 90°), and only
  then explains that **sin-corrected** divides that factor out so flat 1 means random.
- Verified live in both runtimes: Flask + `data/5K_try1` Se–Ga–Se gives 4.00 Se/Ga (4-fold 100.0%)
  at 109.4° ± 3.5°; static Demo (GaTa₄Se₈ 250 K) Se–Ta–Se gives 5.95 Se/Ta (6-fold 95.8%) with the
  split-octahedron distribution. Backend: `TripletsApiTests` (octahedron: exact 12×90° + 3×180°,
  400/403 paths); frontend suite 198 tests incl. the new parity file.
- A 41-agent adversarial review confirmed 17 findings, all fixed: bin-count rounding unified to
  half-up (`_bin_count` — Python's banker's `round()` gave a different **payload shape** than
  `Math.round` at e.g. 8° bins), JS bin assignment now replicates numpy's edge-corrected algorithm,
  dataset switches clear results and epoch-guard in-flight computes, triplet seeding keys on the
  sites payload (not the dataset key racing ahead of it), both app boundaries cap `rmax ≤ 15 Å` and
  `binWidth ≥ 0.05°` and reject half-specified B–C windows (the engine itself stays unrestricted for
  CLI/library use), element case is normalized before the route cache, and assorted UI polish
  (input styling, loading states, debounced helper guides). One documented residual: libm vs V8
  `acos` can differ by 1 ulp, so a *bitwise-ideal* geometry whose cosine sits exactly on a bin edge
  (cos = 0.5 in an undisplaced average configuration) may shift one count between adjacent bins
  across engines — unreachable with real RMC data.

**Bond-angle (triplet) distribution engine + `rmc-triplets` CLI** (2026-08-05) —
`rmc_toolkits/triplets.py` (source of truth, module docstring is the math reference) +
`rmc_toolkits/triplets_cli.py`, following RMCProfile's `triplets_new_bonds_sinth` workflow: pick an
A–B–C triplet with **B the central atom**, bound the A–B and B–C bond lengths (inclusive windows),
and histogram the angle at B. Engine-only by design — offline validation comes before any Flask
route, browser port, or page.

- Neighbours come from a linked-cell search that carries **explicit periodic-image shifts** instead
  of assuming minimum-image: exact for any cell shape (triclinic included) and for boxes smaller
  than the cutoff, where several images of one atom are genuine distinct neighbours. Equivalent
  ends (same element + window) count each unordered bond pair once — an octahedron gives C(6,2)=15
  angles; distinct windows count ordered 1→2/2→3 assignments minus the bond-with-itself
  combinations. All pair formation is vectorized (ragged group cartesian products, no per-atom
  Python loop): the 52 000-atom GaNb₄Se₈ box runs in well under a second per triplet.
- Three curves per histogram: raw `counts`; `density` (per-degree, unit integral); and
  `sin_corrected` — the count fraction over the *exact* isotropic bin fraction
  `(cosθ₁−cosθ₂)/2`, i.e. the "sinth" normalization computed from the bin integral rather than
  `1/sin(θ_center)`, so the 0°/180° bins stay finite and an isotropic configuration reads flat 1.0.
- `rmc-triplets` (console script; `python -m rmc_toolkits.triplets_cli`) takes an `.rmc6f` or a run
  folder, writes a commented CSV (angle, counts, density, sin-corrected) plus an optional PNG and
  raw-angle dump, and prints bond counts, mean lengths, per-central coordination and angle stats.
  Existing outputs are never overwritten without `--force` (the `rmc-autoscale` convention).
- A 14-agent adversarial review (brute-force refutation over ~175 randomized/hand-built
  configurations) confirmed the periodic geometry and counting; its three surviving minor findings
  are fixed: zero-length pairs are never bonds (`rmin=0` + bitwise-coincident atoms used to send a
  0/0 NaN angle past the histogram), the CLI reports unwritable outputs and truncated `.rmc6f`
  headers as clean errors instead of tracebacks, and overwrites now require `--force`.
- `tests/test_triplets.py` (27 tests): exact geometry fixtures (octahedron 12×90°+3×180°), bonds
  through the periodic wall, multi-image bonds in a box smaller than the window, self-image bonds
  of the central element, wrap invariance, **exact agreement with a brute-force all-images
  reference in a skewed triclinic cell**, isotropic sin-correction flatness, normalization and
  validation errors, loader parity, CLI end-to-end. Numpy-on-Accelerate note: large `(P,3)@(3,3)`
  matmuls emit spurious divide/overflow warnings on Apple silicon, so the engine uses `einsum` for
  the displacement→Cartesian map (bit-identical result, no BLAS warning path).
- Physics check on `data/5K_try1` (GaNb₄Se₈): Se–Ga–Se with a 2.1–2.7 Å window peaks at 109.4°
  with 3.99 Se per Ga (GaSe₄ tetrahedra); Se–Nb–Se with 2.2–2.9 Å shows the distorted-octahedron
  splitting (~76/88/104/160°) at 5.78 Se per Nb. Next steps when validation settles: Flask route,
  `workers/` port (parity-tested), and a page.

**Algorithms and math reference** (2026-07-26) — `docs/ALGORITHMS.md` plus seven per-page documents
under `docs/algorithms/`: a code-anchored account of every mathematical operation each page performs
on the user's data, so a reader can audit how a plot, density map, symmetry label, scaled dataset or
direction map was produced. ~16 000 lines.

- The hub carries a page index, a pipeline diagram, a table of *where* each computation runs (Python
  package / Flask route / browser worker), and a consolidated list of every limitation and
  approximation, each linking to the section that explains it. `algorithms/notation.md` reconciles
  the symbol table across sections and lists the symbol *collisions* rather than merging them.
- Two conventions the reference commits to: **the code wins** wherever document and source disagree,
  and *reference-grade* (the float64 Python path) is distinguished from *visualization-grade* (the
  browser ports), with measured cross-engine tolerances quoted — and the gaps named where no such
  test exists (`pcaKde.js` and `orientation.js` are pinned to their own in-language references, not
  to Python; `localKdeWorker.js` has no numerical reference test at all).
- README gains a "Math Under the Hood" section with the signature equation per analysis page — fit
  residual, 2-D Gaussian KDE, ADP tensor with its $\chi^2$ probability ellipsoid, solid-angle
  enhancement — each with the caveat that makes the number honest. QuickStart, REFERENCE and AGENTS
  point at the reference from their own entry points.

**Rwp reports "unavailable" instead of a fake perfect fit** (2026-07-25) — `rwp()` returned `0`
whenever its denominator was falsy. For an observed column that is entirely NaN (a fully masked
dataset) the denominator is NaN, `!NaN` is `true`, and the dashboard chip showed `Rwp 0.000` — read
as a perfect fit over data that isn't there. A genuinely zero denominator collapsed to the same `0`.

- Both implementations now sum only the points where the observed *and* fitted values are finite,
  and return `None` / `null` for the two degenerate cases (no finite pair, zero denominator).
  `rmc_toolkits/parsers.py` stays the source of truth; `browserData.js` is the static-mode port —
  keep them in sync (see AGENTS.md). `PlotResult.metrics` is now `dict[str, float | None]`, so
  `/api/plot/metadata` serializes the sentinel as JSON `null`, matching static mode.
- `Dashboard.jsx` renders the sentinel as `Rwp —` via a single shared `renderRwpChip` helper (the
  two chip sites had duplicated the formatting). The assistant's run context already skipped
  non-finite Rwp values, so it drops the dataset's `rwp` key as before.
- Tests: all-NaN observed, partially-NaN columns on either side, and a zero-denominator column, in
  both `tests/test_parsers.py` and `src/__tests__/browserData.test.js`.

**PCA Ellipsoid: principal axes in the crystallographic frame** (2026-07-25) — the Displacement
statistics panel gains a fourth column, **Crystal orientation**: for every principal axis, its angle
to the unit-cell edges a, b, c (the closest one shaded, same accent wash the covariance diagonal
uses) and the crystallographic direction `[u v w]` it runs along. `src/pcaCrystalFrame.js` had
computed all of this since v0.4.0 but nothing imported it — only `unitCellVectors` was wired up, so
the app knew which lattice direction each displacement axis pointed along and never said so.

- Angles are between Cartesian directions (both frames already share one Cartesian basis), so they
  are well defined for any cell, oblique or not. `[u v w] = M⁻ᵀ·axis`, normalised so the largest
  component is 1, is a **direct-lattice** direction — in an oblique cell it is not normal to the
  like-indexed (h k l) plane. The panel's InfoBadge says both, since the distinction is easy to
  misread off a table of numbers.
- An eigenvector's sign is arbitrary (v and −v are one axis), so the new
  `crystalOrientationRows` reports the member of the ± pair whose closest crystal axis is the
  *acute* one — angles and `[u v w]` flip together, keeping each row internally consistent, and
  "dominant" then literally means smallest angle. Unit-tested in `pcaCrystalFrame.test.js`
  (already-acute rows untouched, a 170°-from-a axis reported as 10° with the direction negated, no
  obtuse dominant angle for an oblique cell, singular cell → null).
- The column appears only when the file carries lattice metadata; the statistics grid then switches
  to four tracks (`.has-crystal`, sized to avoid inner scrollbars down to a 1280-wide screen) and
  still stacks to one column below 980 px.
- Also: the header wordmark centers "WORKBENCH" under "RMCProfile" instead of leaving it flush-left,
  so the two lines read as one lockup.

**Displacement-orientation engine (hex-binned sphere)** (2026-07-24) — new engine for the
distribution of displacement *directions* per site: `u = Δr/|Δr|` binned in solid angle on
a Goldberg polyhedron (the dual of a frequency-ν geodesic icosahedron — hexagons plus the
12 pentagons any hexagonal tiling of a sphere must carry; 10ν²+2 cells), each cell's count
divided by its exactly-computed solid angle. Complements the ellipsoid view: discrete
hop-sites and antipodal asymmetry (odd anharmonicity) are visible here and invisible in U.

- `rmc_toolkits/orientation.py` (source of truth) + `workers/orientation.js` (static-mode
  port, keep in sync): tiling construction with derived (not hard-coded) icosahedron faces
  and combinatorial-key merging, O(1) direction→cell lookup (gnomonic seed + greedy walk,
  brute-force-verified), exact centrosymmetric antipode map.
- Histogram options: auto resolution via `recommended_frequency` (over-binning guard,
  ~12 pts/cell), weights `count | amplitude | amplitude2` (the latter's angular map is the
  decomposition of ⟨u²⟩), amplitude cutoffs (absolute + quantile — near-zero |Δr| has
  noise for a direction), mass-conserving neighbour smoothing, `cartesian | pca` frame
  (PCA axes shared with the ellipsoid engine via the same canonicalization). The map is
  deliberately never antipodally folded: the +u/−u imbalance is the signal.
- Outputs: `enhancement = 4π·density` (1 = isotropic, reads as "N× more likely than
  chance"), per-cell Poisson `zScore` vs the isotropic null (smoothing never feeds the
  z), `antipodalAsymmetry` = Σ_pairs|n(u)−n(−u)|/N with its Poisson noise floor
  `antipodalAsymmetryNull` (inversion asymmetry U is blind to), orientation tensor
  ⟨u uᵀ⟩ + eigenvalues + `orientationAnisotropy` (3λ₁−1), cell polygons for rendering.
- Wiring: `/api/pca/orientation` (Flask), `{kind: 'orientation'}` in `pcaKdeWorker.js`,
  package exports. Tests: `tests/test_orientation.py` and
  `workers/__tests__/orientation.test.js` — tiling invariants, brute-force assignment
  parity, isotropic flatness, lobe recovery, one-sided asymmetry detection, weight/frame
  behavior, API + worker routing.
- **Renamed the nav tab "Orientation" → "Displacement Directions"** (2026-07-24, same
  branch): clearer, and the direct direction-counterpart to "PCA Ellipsoid" (amplitude /
  shape). The internal route key stays `orientation`. Also: the Axis-views hover now
  darkens with a real ~38% scrim (a weak dark overlay over the light panel just greyed
  toward white, and the global `button:hover` had out-specified it), and the a/b/c
  labels are larger.
- **Orientation controls: contrast knob + defaults** (2026-07-24, same branch): new
  **Contrast** slider (0.5–3×) applies a symmetric color gain about the isotropic level
  (`colorCoordinate` in `orientationSphere.js`: `t = clamp(pivot + contrast·(v/vmax −
  pivot))`, `pivot = 1/vmax`), so faint lobes and depletions stand out; contrast = 1 is
  exactly the old linear mapping (backward compatible), and `colorbarGradient` paints
  each stop through the identical transfer so the bar and the sphere always agree
  (unit-tested). Defaults changed to **ν = 10** (was Auto) and **2× smoothing** (was 0)
  for a legible out-of-the-box map, and the smoothing slider now reaches **12×** (was 4).
  The redundant "N cells (ν=…) / N vectors" summary line is removed (both are already in
  the Resolution control and the panel header). Axis-views hover is now a translucent
  accent tint + ring instead of a near-white wash.
- **Orientation page — three-panel layout + PCA-page design parity** (2026-07-24, same
  branch): the page now matches the PCA Ellipsoid page's conventions — all options live
  in the top controls bar, and the header actions (Crystal|PCA frame toggle, Reset view,
  Save) sit in each panel's own header. Three equal-height panels in one grid at
  **3 : 6.5 : 6.5** — Axis views (the fixed-angle mini column, now its own panel) :
  sphere : site picker; `OrientationView` renders `display:contents` so its two panels
  drop straight into the grid. Panel height is viewport-clamped
  (`calc(100vh − 17.75rem)`) so the three fill a 16:9 screen down to the footer with zero
  scroll (verified 1600×900). The mini views zoom out (camera pushed to 5.6×) so a
  fully-inflated relief surface never clips the pane edges, and the main sphere's default
  view is pulled back too. The axis rods now show **only the selected coordinate frame**
  (a/b/c in Crystal, PC1/2/3 in PCA) — the mini views, main sphere, and legend all show
  one system at a time instead of overlaying both.
- **Orientation page 16:9 layout + shared site picker** (2026-07-24, same branch): the
  clickable Site-ellipsoids unit-cell picker is extracted from PcaKdePage into
  `SiteStructurePanel.jsx` (axis palettes/builders into `sceneAxes.js`) and now sits on
  the Orientation page too — sphere panel wide left, picker right
  (`.orient-page-layout` grid; the panels' named grid-areas from the PCA layout are
  neutralized there). The fixed-angle mini views moved from a bottom row to a **column
  left of the sphere** (scissor panes stacked vertically; WebGL y counts from the
  bottom). The sphere row height is viewport-clamped (`calc(100vh − 30.5rem)`, self-sized
  `flex: none` so a zero-height flex ancestor can't collapse it) — on a 16:9 monitor the
  whole page (site bar, controls, sphere + minis + picker, readouts, footer) fits with
  zero scroll, verified at 1600×900.
- **Orientation promoted to its own workspace page** (2026-07-24, same branch): the
  histogram is not a PCA product, so the Density | Orientation tab coupling is gone —
  `OrientationPage.jsx` is a top-level nav page with its own site picker, and the PCA
  Ellipsoid page is back to its pre-tab layout. The shared plumbing (text loading,
  worker/API routing, sites table, selected site) moved to the `useSiteCloud` hook with
  one app-lifetime worker, so both pages share the parse cache. The sphere view gains
  **three fixed-angle mini views** (down a/b/c, or PC1/2/3 in the PCA frame): one extra
  scissored renderer re-rendered only on scene changes (static cameras, no per-frame
  cost); clicking a mini snaps the main camera to that axis view.
- **Amplitude relief** (2026-07-24, same branch): the engines report
  `cellMeanAmplitude` — the mean |Δr| of each cell's movers (smoothing applies to the
  numerator/denominator sums, not the ratio; empty cells report 0) — and the sphere gains
  a **Relief** slider that bulges each cell radially by its mean |Δr| relative to the
  site average. Color (how often) and shape (how far) carry independent information.
  Shared polygon vertices average their cells' factors so the relief surface is
  crack-free (unit-tested); cell borders follow the surface; the hover readout adds
  ⟨|Δr|⟩.
- **Adversarial review pass** (multi-agent, 2026-07-24) — confirmed findings fixed:
  auto-resolution rounding divergence (Python `round()` is half-to-even, JS used
  `Math.round` half-up — at exactly 774 surviving points the two engines picked
  different tilings; JS now rounds half-to-even and both suites pin
  `recommendedFrequency(774) == 2`); stale hover index could crash the tooltip when a
  resolution change shrank the tiling (hover now resets with each result + bounds
  check); OrientationView unmount leaked its geometries and a WebGL context per tab
  switch (groups disposed + `forceContextLoss` in the scene cleanup); the
  density-integrates-to-one tests were circular (areas cancel by construction) and are
  replaced by a Monte-Carlo cross-check that assignment Voronoi fractions × 4π
  reproduce the analytic polygon areas. Greedy-walk assignment verified exact against
  brute force at ν = 24/40/64 (0/20000 mismatches, all centers self-assign).
- **UI: Orientation tab** in the PCA Ellipsoid main panel (`OrientationView.jsx` +
  `orientationSphere.js` pure helpers, unit-tested), sharing the page's site picker.
  Density | Orientation tabs in the panel header (the density canvas hides via CSS, never
  unmounts — its once-mounted scene effect survives tab switches). Flat-shaded Goldberg
  cells colored by enhancement, cell borders, PC1/2/3 + a/b/c axis rods, Crystal ↔ PCA
  frame toggle, per-cell hover readout (direction, ×isotropic, count, z), colorbar with a
  1× tick, and a summary strip: cells/ν, used vectors, peak (+z), orientation anisotropy,
  ± asymmetry vs its Poisson floor (flagged red when > 3× floor), map significance.
  Verified live in both runtimes (Flask `data/RMC`, browser Demo run).

**Auto StoG tab hidden** (2026-07-23) — the page is not behaving correctly yet, so it is
gated out of the UI while the work continues offline. `App.jsx` gains a `SHOW_AUTO_STOG`
constant (currently `false`) guarding the nav button and the page mount; with it off the
page is unreachable, since `activePage` defaults to `dashboard` and that button was the
only route to it. Nothing was removed — `AutoStogPage.jsx`, `workers/autoScale*.js`,
`rmc_toolkits.scaling`, the `rmc-autoscale` CLI, and `/api/scaling/*` are untouched, and
their tests still run. Flip the constant to restore the tab.

**Auto StoG independence + ρ₀ self-consistency** (2026-07-18) — the page is now true
pre-processing, decoupled from the run folder, and the number density can come from the
data itself.

- **ρ₀ self-consistency** (`estimate_rho0` in `scaling.py`, JS port in `autoScale.js`,
  CLI `--estimate-rho0`): the density-limit amplitude depends on ρ₀ (C2 rows scale with
  the density line −4πρ₀r) — degenerate with the scale from low-r alone — while the Q→0
  Faber-Ziman amplitude needs only ⟨b²⟩ + the measured level. ρ₀ is recovered as the root
  of `a_fz/a_density(ρ₀) = 1` by fixed-point iteration (`ρ ← ρ·concordance`, 2–4
  `autoscale` passes). Requires a composition; `extrapolated` flags Qmin beyond the FZ fit
  width. Validated: synthetic truth 0.05 recovered within 2% from seeds 0.02–0.2; FeCoSn
  x-ray 199 K gives 0.0600 vs the hand 0.057329 (<5%) from any seed. ρ₀ resolution order
  is now: value → `NUMBER_DENSITY ::` header → mass density + composition →
  self-consistent estimate (auto-run when ρ₀ is left empty; explicit "Estimate ρ₀" button
  too) — so composition + Q window are the only required inputs.
- **Auto StoG page decoupled from the run folder** (`AutoStogPage.jsx`): pre-processing
  now has page-local file handling — an upload/dropzone (S(Q) ± stog.inp, multi-file;
  stog output files filtered out) replaces the shared run-folder listing, and the page is
  fully client-side in BOTH runtimes (worker engine + zip export; the Flask
  `/api/scaling/*` endpoints stay for API/CLI use but the page no longer calls them).
  Q window prefills from the data's finite extent; `App.jsx` passes no props.
- **Parameters grouped with descriptions**: the Advanced panel is five fieldsets, each
  with a one-line purpose and per-field tooltips — *Amplitude & offset* (High-Q mode,
  amplitude criterion, Robust/σ/Despike), *Coefficients* (⟨b⟩², ⟨b²⟩ x-ray overrides),
  *Transform* (r-cut, grid, Lorch, low-Q corr.), *Low-r region* (r₀, fit window,
  enforcement), *Fixed scaling* (a, b).
- Tests: `estimate_rho0` synthetic + x-ray + CLI coverage (Python), fixture-backed JS
  parity (tolerance bounded by the rtol stopping rule, not single-pass precision);
  browser flow verified end-to-end in static mode (upload → estimate 0.060664 →
  auto-scale a = 1.2511, concordance 1.00 → 9-file zip).
- **Correctness review of the page** (multi-agent sweep of every description + flow +
  export against the engine; 12 findings, 8 confirmed after adversarial verification):
  manual runs now recover r₀ by detection so a checked "Enforce low-r" with an empty
  cutoff is honored (was a silent no-op and a CLI-parity break); a fixed-(a, b) run no
  longer fails the unused FZ-amplitude validation; "Run fixed" gains the missing
  `estimating` concurrency guard; three tooltips corrected to the code's real resolution
  orders (ρ₀ chain includes stog.inp above the data header; r₀ is header → stog.inp peak
  window → detection; x-ray needs ⟨b⟩² = 1 *and* ⟨b²⟩); export title lines use %g-style
  numbers (the `(browser)` marker vs the CLI's version token is intentional); provenance
  JSON now records `stogInpReference` (the loaded stog.inp's hand scaling), the fit
  `history`, and `c1ModeEffective`, matching the CLI payload.
- **Discordant-data guardrails** (found live on the Mn₃Sn neutron runs, whose low-Q
  hole breaks the density-limit amplitude): `estimate_rho0` stops on non-physical
  amplitudes instead of iterating to a bound, and an unconverged estimate is *refused*
  everywhere (CLI error, worker error, page banner) with guidance to set ρ₀ explicitly
  and use the FZ amplitude — previously the page silently fit with a garbage density
  (a = −1.02 on PG3_54139). The neutron-vs-x-ray state is now always visible: a chip
  shows the coefficients in effect and their source, warns when Advanced overrides
  shadow a typed composition, and a "Reset params" button clears a previous sample's
  settings. Folder drops report cleanly instead of failing silently. Engine + CLI
  guard tests (Mn₃Sn-fixture-gated).

**Composition-first Auto StoG** — procedure clarified end-to-end
([SCALING_PROCEDURE.md](SCALING_PROCEDURE.md)); the user now provides only the chemical
composition + [Qmin, Qmax] (+ a density when the data header has none); validated on three
new complete Fortran neutron runs (POWGEN Mn₃Sn).

- **The procedure document** (`docs/SCALING_PROCEDURE.md`): the definitive recipe for
  absolute-scale S(Q)/G(r) — inputs and defaults tables, the Keen-convention function map,
  the 7-step pipeline, how to read the one-sided density flag and the concordance trust
  metric, and material-class guidance. Research base: pystog 0.6.7 source (operation order
  confirmed read → merge → manual scale → transform → filter → Lorch → Keen conversions;
  no auto-scaler, low-r minimizer an explicit TODO) and ADDIE (density conversions,
  periodictable molar masses, `<b_coh>^2` hand-off).
- **Composition-aware omitted-low-Q correction**: the analytic [0, Qmin] correction now
  extrapolates S(Q) to the composition-derived Keen Eq. 21 target
  ``S(0) = 1 − ⟨b²⟩/⟨b⟩²`` instead of pystog's solid-state 0 — algebraically a
  one-line change (``const' = (1 − s0)·const + s0·coef``, still affine in (a, b)). For
  negative-b compositions this is O(1): Mn₃Sn has ⟨b²⟩/⟨b⟩² = 13.06 → S(0) = −12.06, and
  the composition-aware target cuts the low-r residual ~40% on the PG3 runs.
  ``ScalingConfig.s0_target`` (auto from ``b_sq_avg``; pin 0 for classic parity).
- **Data-derived first-shell r₀** (`detect_first_peak_onset` + a second refinement pass in
  `autoscale`): the dominant |g| feature's left flank (35% of peak) — |g| because
  negative-b totals can have an *inverted* first shell; peak-relative because sub-r₀
  ripples scale with the amplitude. Detected 2.73–2.77 Å on the Mn₃Sn runs and 2.53 Å on
  FeCoSn (hand-chosen classic cutoffs: 2.40–2.68). Sets the low-r fit window without a
  structural prior and becomes the default low-r **enforcement cutoff** in data mode
  (CLI + API + page; explicit `--no-enforce`/`enforce:false` still opts out);
  `diagnostics_summary` reports `r0_detected`/`window_refined` and the *effective* fit
  window from provenance.
- **ADDIE-style density toolkit** (`scattering.py` + JS port): `ATOMIC_MASS_U`
  (periodictable/CIAAW values for the full element set), `molar_mass`,
  `number_density_from_mass_density` / inverse (ρ₀ = ρ_m·N_A/10²⁴·n/M). CLI
  `--mass-density`, API `massDensity`, page "or ρ g/cm³" field. Mn₃Sn check:
  0.063049 atoms/Å³ ↔ 7.421 g/cm³.
- **Mn₃Sn neutron validation** (`data/stog_tests/{stog,stog_300K,stog_500K}`, local):
  three complete Fortran runs; parameters recovered from the outputs to ~1e−12
  (hand scalings ×2.5/×2.05/×10 — all satisfying b = 1 − a, i.e. the level-anchored
  decomposition the sweep formalizes; ⟨b⟩² 0.015407 = `faber_ziman("Mn3Sn")` exactly).
  New `Mn3SnNeutronTests` (+ detection/s0 unit tests): manual filter-stage parity ≤2e−3
  rms per run, composition/density round trips, first-shell detection, and the headline
  physics finding — the density limit is degenerate on this material (one-sided flag
  False, hand values mutually inconsistent by 5×) while the **FZ amplitude injects the
  composition's S(0) and lands at O(10)** consistently; the procedure doc encodes that
  decision rule.
- **Page**: primary bar is now Composition · Qmin · Qmax · ρ₀/mass-density with a live
  ⟨b⟩²/⟨b²⟩/S(0) chip; ⟨b⟩²/⟨b²⟩ overrides moved to Advanced (x-ray note kept); a
  collapsible "How Auto StoG works" explainer; a First-shell-r₀ readout card; and
  **full-range G_K(r) and D(r)** (no 8 Å slice — theory guide lines stay confined to the
  low-r region; the G_K default y-zoom keeps the −⟨b⟩² level readable, box-zoom for
  detail). JS engine/worker fully synced (s0-aware correction, detector, two-pass,
  converters) with new Python-golden parity tests (detection case + Mn₃Sn converters).
  Verified live on the 300 K run: composition-only inputs → a = 1.309, 7 iterations,
  r₀ 2.750 Å detected + enforced, density flag red, concordance 7.95 flagged.

Auto StoG page redesign + **Phase 5: the static-mode engine** — the hosted app now runs
Auto StoG entirely in the browser.

- **Static-mode Auto StoG (plan Phase 5)**: `src/workers/autoScale.js` is a straight JS
  port of the Python engine — transforms (trapezoid sine FT, Lorch, omitted-low-Q
  correction, Fourier filter, first-peak zeroing), the level sweep (prefix-sum OLS with
  numpy `.astype(int)` edge parity), the Huber-IRLS closed-form fit, sweep/joint/FZ
  amplitude modes, despiking, diagnostics, the stog parsers (`stog.inp`, xy files, `::`
  headers, Fortran-style writer), and the Faber-Ziman calculator (Sears table + formula
  parser). `autoScaleWorker.js` runs it off-thread with transferable buffers. **Parity is
  tested, not assumed**: `tests/generate_autoscale_fixture.py` freezes Python golden
  numbers into `src/__tests__/fixtures/`, and `autoScale.test.js` (11 vitest tests)
  asserts the JS fit matches to 1e−6 relative (level sweep to 1e−9, manual pipeline
  samples to 1e−9). Verified live: the in-browser fit on the FeCoSn 199 K run reproduces
  the backend exactly (a = 1.1839, b = −0.19804, 3 iterations, level 1.0120 ± 0.017).
  Static export packs the classic 9-file family into a zip (client-side `writeStogXy`,
  provenance JSON included); the Auto StoG tab now shows in both runtimes and
  `browserData.isSupportedFile` admits `.inp`.
- **Page redesigned on the app's design language, laid out for 16:9** (was hardcoded-hex
  cards in a narrow sidebar): a PcaKdePage-style horizontal controls bar (SOURCE picker +
  micro-labeled parameter fields + Auto-scale + Advanced pill), an expandable advanced
  bar (windows/grid, sweep-vs-joint, density-vs-FZ, toggle pills, fixed-(a, b) expert
  run), a stat-card readout strip (correction vs hand values, convergence with the
  per-iteration a-trajectory, high-Q level ± uncertainty, fit quality, density-limit
  verdict, concordance), then a full-width S(Q) card over side-by-side G_K(r)/D(r)
  cards — all on `index.css` tokens (borders, shadows, accent, pills, tabular numerals).
- **InteractivePlot gains guide-line support** (additive; Dashboard untouched): series
  with `role: 'guide'` render dashed/muted outside the palette rotation and are skipped
  by hover snapping (dashed legend swatches); `defaultHidden` series start muted;
  `initialYDomain` sets the un-zoomed default view (used to keep the G_K low-r level
  readable instead of the first peak); swapping `plotData` now resets zoom/hover/hidden
  state (previously stale zoom survived a re-fit).
- **Plot content**: measured-unscaled S(Q) ships default-hidden (one legend click away),
  the S → 1 asymptote, the measured level (drawn only over its admissible window), and
  the S(0) Faber-Ziman target marker join S(Q); theory lines −⟨b⟩² and −4πρ₀⟨b⟩²r anchor
  the G_K/D(r) cards.
- **Form & workflow**: session persistence (source + all settings survive a reload via
  sessionStorage), the Formula field is labeled *(neutron)* with an x-ray tooltip, an
  inline guard disables Auto-scale when FZ mode lacks ⟨b²⟩, and micro-labels no longer
  pass through `text-transform: uppercase` (which had turned ρ₀ into a capital-rho
  P-lookalike).

Auto StoG Phases 3–4 — the Flask scaling API and the **Auto StoG** page.

- **`/api/scaling/preview` + `/api/scaling/run`** (`web_app/backend/app.py`): POST endpoints
  driving the Phase-1 engine server-side. `preview` resolves a classic `stog.inp` or a bare
  data file (+ overrides; `.dat`-header ρ0/r0 and `formula` defaults mirror the CLI), runs the
  auto-fit (or a fixed manual scaling) behind a per-(path, mtime, config) LRU cache, and
  returns the full series (raw/scaled/filtered S(Q), G_K, D(r), enforced variants), theory
  guide values, diagnostics, and provenance; `{"inspect": true}` is the cheap no-compute form
  used to pre-fill the page. `run` writes the classic output family through the *same* writer
  as the CLI (`ft.dat` included) with the no-clobber guard mapped to HTTP 409, outputs
  restricted to the configured data roots. 5 new backend tests (synthetic run under
  `results/`; 20-test backend module green).
- **The Auto StoG tab** (`AutoStogPage.jsx` + CSS): first in the tab row (Dashboard remains
  the default page), Flask mode only — static mode hides the tab and the page shows a
  pointer to the local app. Automation-first per the plan: pick a source file (stog.inp
  pre-fills everything; data files pre-fill from the `.dat` header), one **Auto-scale**
  button, then a diagnostics readout (a/b beside the stog.inp hand values, convergence +
  level ± uncertainty, low-r rms, one-sided density-limit verdict, amplitude concordance
  with a "check ρ₀" hint on discord) over three InteractivePlot charts with theory
  guide-lines — S(Q) (asymptote + measured level), G_K(r) (−⟨b⟩²), D(r) (−4πρ₀⟨b⟩²r) —
  the latter two zoomed to the low-r region, enforced curves overlaid when enforcement is
  on. Advanced panel: r-windows, Lorch/despike/robust/low-Q-correction/σ toggles, sweep vs
  joint architecture, density vs FZ amplitude criterion, and fixed-(a, b) expert runs.
  Export card writes the RMCProfile-ready family (default `autoscale/` beside the input,
  Force required to overwrite, 409 surfaced as a hint). Verified end-to-end against the
  FeCoSn 199 K run in the live app: page auto-fit reproduces the CLI exactly
  (a = 1.1839, b = −0.1980, level 1.0120 ± 0.017, concordance 1.06).
- File browser (`/api/files`) now also lists `*.inp`, `*.sq`, and `*.dat` so scaling
  sources are pickable.
- Docs-on-completion pass: ROADMAP Phase 7 "Deferred preprocessing" → **active Auto StoG**
  feature; AGENTS.md architecture map gains the four scaling modules, the endpoints, and
  the page. Remaining stretch: plan Phase 5 (static-mode Web-Worker port).

Auto StoG Phases 1–2 — automatic total-scattering data scaling engine (`rmc_toolkits.scaling` +
`rmc_toolkits.transforms`) and the `rmc-autoscale` CLI, plus the app rename to
**RMCProfile Workbench**.

- **Faber-Ziman amplitude mode** (idea: Tsung-Han Yang): `ScalingConfig(
  amplitude_criterion="fz")` / CLI `--amplitude fz` implements the "subtract the measured
  high-Q level, scale Q→0 onto S(0) = 1 − ⟨b²⟩/⟨b⟩², shift the level back to 1" procedure —
  closed form on top of the level sweep, no self-consistent loop (the criterion is
  filter-independent), requires ⟨b²⟩. Because it never touches ρ0, it is the natural
  cross-check for the density-limit amplitude: on the FeCoSn validation data a ±10% ρ0
  error moves the density amplitude ~1:1 while the fz amplitude is bit-identical, so the
  concordance diagnostic turns a wrong `NUMBER_DENSITY ::` into a measurable discord. In
  fz mode `diagnostics_summary` suppresses the (vacuous) self-concordance and the density
  residuals act as the independent check.
- **FeCoSn 199 K validation + robustness campaign** (`data/stog_tests/199K`, local; script
  `data/stog_tests/robustness_199K.py` + results JSON): a third complete classic-Fortran
  run now validates the stack — manual parity `scale.fq` 6.8e−14 max|Δ|, `ft.dat`/
  `scale_ft.sq` 2.1e−5 rms, enforced `scale_ft_rmc.gr` 9.6e−5 rms. 63-case perturbation
  study: affine pre-corruption invariance exact (≤2.6e−11%); recovered scale stable ±3%
  over Qmax ∈ [18, 28] and Qmin ∈ [0.5, 1.6]; noise graceful (±0.3% at σ = 0.005, ±2.6%
  at σ = 0.05); spikes −19% unflagged→flagged, `--despike` restores to +4.3%; every
  catastrophic case (rolloff Qmax, starved Qmin, spikes) raised
  `density_limit_satisfied=False`. All three auto criteria (sweep+density, joint, fz)
  land 8–13% above the colleague's hand scaling while improving the honest low-r
  residual (auto 0.110 vs hand 0.165, −33%), with density/fz concordance 4.6%. The
  x-ray fixture tests now select the first available FeCoSn run (`100K` or `199K` — same
  stog parameterization), and a CLI-level x-ray parity test joins the suite. 139 tests.
- **Auto StoG Phase 2 — the `rmc-autoscale` CLI** (`rmc_toolkits/scaling_cli.py`; console
  entry via `[project.scripts]`, module form `python -m rmc_toolkits.scaling_cli`): drop-in
  replacement for an interactive classic-stog session. Reads a classic `stog.inp` — or
  `--data FILE --qmin --qmax` with `--formula`-computed coefficients (scattering.py) and
  ρ0/r0 pre-filled from the `.dat` `NUMBER_DENSITY ::`/`MINIMUM_DISTANCES ::` header — and
  auto-fits (a, b) by default; `--manual` reruns the stog.inp hand scaling and
  `--scale`/`--offset` fix them explicitly. Writes the seven classic output files (scaled
  S(Q), unfiltered g−1, filtered S(Q), filtered g−1 + D(r) companion column, and the RMC
  FK/GK/D(r)) plus a provenance JSON carrying the full configuration, fit history, and
  `diagnostics_summary`. Safety per plan §5: outputs default into an `autoscale/` directory
  beside the input and nothing is overwritten without `--force`, so the tool can never
  silently clobber the real STOG outputs a `stog.inp` sits beside. The RMC files get the
  exact Fortran `first_peak_zero` enforcement by default in stog.inp mode (`--no-enforce`
  opts out); the honest pre-enforcement low-r rms is always printed. `write_stog_xy` gains
  an optional third column for the `scale_ft.gr` layout. The classic fixed-name `ft.dat`
  Fourier-filter correction is written too (data-grid), so the output family matches a
  stog/pystog session file-for-file. Tests: `tests/test_scaling_cli.py`
  (9: synthetic auto/manual/data-mode end-to-end, no-clobber guard, error surfaces,
  module-entry smoke, skip-if-absent Fortran parity). Full suite: 135 tests green.

- **New `rmc_toolkits/transforms.py`**: Keen-2001-convention conversions (S(Q) ⟷ F(Q) ⟷ F_K(Q);
  g(r) ⟷ G_PDF(r) ⟷ G_K(r) ⟷ D(r)), the trapezoid sine-FT pair, Lorch window, the analytic
  omitted-low-Q correction (affine-basis form), the classic stog/pystog Fourier filter, and the
  stog low-r enforcement stage. Discretization validated against pystog 0.6.7 and a complete
  classic Fortran stog run (`data/stog_tests/stog_59438`, local-only): filter correction and
  filtered S(Q) agree to ~6e-4 rms; enforced RMC outputs match the Fortran files to 1e-9.
- **New `rmc_toolkits/scaling.py`**: the auto-scaler. Affine correction `S_corr = a·S + b`
  (multiply convention; classic stog's `yoffset/yscale` map via `a = 1/yscale`), fitted by a
  closed-form linear least-squares against the high-Q asymptote (Keen Eq. 21) and the low-r
  density limit (Eqs. 15/29 in g-space), inside a self-consistent loop with the Fourier filter.
  Recovers known (a, b) on synthetic data to ~0.3% (the omitted-low-Q correction, on by
  default, is what makes that possible — 8% bias without it). `diagnostics_summary` reports the
  honest pre-enforcement residuals and flags datasets whose missing low-Q information makes the
  absolute scale unrecoverable from self-consistency (the 59438 example is such a case: its
  filtered outputs violate the Krogh-Moe sum rule ~26× even at the expert's hand scaling).
- **Parsers**: `StogInput`/`read_stog_inp` (classic 23-line stog.inp, with explicit
  `NotImplementedError` on unexercised variants), `read_stog_xy` (count-header/NaN-tolerant),
  `read_dat_header` (`TITLE ::` / `NUMBER_DENSITY ::` / `MINIMUM_DISTANCES ::`), and
  `write_stog_xy`.
- **Tests**: `tests/test_transforms.py` + `tests/test_scaling.py` (30 tests): synthetic
  round-trips and known-scale recovery always run; Fortran-run parity tests skip cleanly when
  the local example is absent. Full suite: 95 tests green.
- **New `rmc_toolkits/scattering.py` — Faber-Ziman coefficient calculator**: bound coherent
  neutron scattering lengths for 89 natural elements (NIST NCNR / Sears 1992; real part for
  the complex-b absorbers B, Cd, Dy, Eu, Gd, In, Sm), a chemical-formula parser (decimals,
  parentheses: `"Sr0.5Ba0.5TiO3"`, `"Al2(SO4)3"`), and `faber_ziman()` returning ⟨b⟩² (the
  stog "Faber-Ziman coefficient") and ⟨b²⟩ in both barns and fm² — the ecosystem mixes units
  (pystog's argon example is fm²; classic stog inputs are barns). Cross-validated against
  pystog's argon config (3.644 fm², exact). Per-element overrides support isotopic samples;
  null-matrix compositions (⟨b⟩ ≈ 0) are rejected with a clear error.
- **App rename**: "RMCProfile Run Monitor" → **"RMCProfile Workbench"** (header, tab title,
  READMEs), reflecting the multi-tool scope ahead of the Auto StoG page.
- **Adversarial review hardening** (15-agent verified review; 11 confirmed findings fixed):
  Q ≤ 0 grid rows are cropped and `fourier_filter` rejects non-positive grids (NaN poisoning);
  the omitted-low-Q correction returns exactly zero when data start at Q = 0 (pystog parity —
  no double counting); the Lorch-branch removable singularity at r = π/Qmax gets its analytic
  limit; `np.trapezoid` shim for NumPy < 2.0; `nr`/`rmax`/`yscale` validation;
  `read_stog_xy` picks the modal column count (numeric headers can't eat the data);
  `fk_to_sq` exported. Most important: the diagnostics flag is now the one-sided
  `density_limit_satisfied` — verification *demonstrated numerically* that a smooth missing-
  low-Q deficiency is silently absorbed into a ~9–21% biased scale with all residuals clean,
  so **False proves the absolute scale is unrecoverable, but True does not certify it**.
- **Level sweep — a criterion-driven answer to "what Q is high enough?"** (idea:
  Tsung-Han Yang). `level_sweep()` searches every candidate high-Q window (both edges swept,
  O(1) per-window fits via prefix sums); a window is *admissible* iff its slope is
  statistically zero given its own fit noise (no hand-set tolerance), the minimum-variance
  admissible window wins, and the level spread across all admissible windows is the honest
  level uncertainty. End artifacts exclude themselves — the criterion independently
  rediscovered both experts' hand cuts (FeCoSn: 24.5 vs hand 26 with the rolloff onset
  caught earlier; PG3: 28.8 vs hand 28). `autoscale` now defaults to the **sweep-anchored
  architecture** (`c1_mode="sweep"`): offset tied to the measured level (`b = 1 − a·level`),
  leaving the density limit a single amplitude dof — the "shift by the level, then scale"
  decomposition, which removes the 2-dof level/amplitude trade-off pathologies and converges
  ~4× faster. No flat window → `asymptote_found=False` and automatic fallback to the joint
  fit. Reported in provenance and `diagnostics_summary`.
- **Dual amplitude criteria with a concordance diagnostic**: alongside the density-limit
  amplitude, `amplitude_from_fz_limit()` independently estimates the scale from the Q→0
  Faber-Ziman limit (Keen Eq. 21: `S(0) = 1 − ⟨b²⟩/⟨b⟩²`, robust low-Q extrapolation;
  requires `b_sq_avg`). `diagnostics_summary` reports both, their ratio, and an
  `amplitudes_concordant` flag — agreement is evidence the absolute scale is trustworthy;
  disagreement *quantifies* how much the data cannot decide it (FeCoSn: 12% discord).
- **Robust high-Q level fitting**: Huber IRLS re-weighting of the joint fit (default on;
  per-block MAD scaling so C1 and C2 are each protected against isolated outliers), optional
  per-point `sigma` 1/σ-weighting of the high-Q rows, an experimental `c1_slope_nuisance`
  term absorbing linear tail drift in the level estimate, and opt-in rolling-median
  `despike` for detector-glitch contamination — measured to restore clean 0.3% recovery
  under tail spikes that otherwise ring through the transform into the low-r window (a
  channel row re-weighting cannot reject: ~80% scale error without despiking). Despike stays
  OFF by default because it also flags real Bragg maxima on crystalline data (12% of points
  on the 59438 benchmark); the dropped count is reported in provenance (`n_despiked`).
- **Second validation dataset — FeCoSn 100 K x-ray run** (`data/stog_tests/100K`, local):
  exercises the normalized-S(Q) conventions (⟨b⟩² = 1, flat −1 level, `1.0 0 0`
  enforcement). Fortran parity to 7.9e−6 rms through the filter stage; the auto-scaler lands
  within 8% of the expert's hand values and improves the low-r residual by 26% — the
  density limit is satisfiable here (Qmin = 0.5 Å⁻¹), unlike the neutron 59438 case, and the
  one-sided diagnostic correctly reports True.
- **Exact Fortran final-step semantics** (from `stog_new3.f90`, located during verification):
  `first_peak_zero()` implements the real ripple-removal rule — zero g(r) where r ≤ cutoff
  *and* outside the first-peak window [rmin, rmax] — which degenerates to the flat −⟨b⟩²
  replacement for the validation example's parameters. The mysterious second stog.inp
  "yoffset" is the Fortran's global "Add values" knob (`y·(1+vadd)+vadd`); still rejected as
  unsupported when nonzero.
- Plan: [STOG_SCALING_PLAN.md](STOG_SCALING_PLAN.md) (build phases, verified math spec,
  validation results).

## v0.4.0 — 2026-07-16

PCA Ellipsoid: PC ⟷ crystal reference-frame switch, and Site-ellipsoids crystal axes + reset/save.

- **The main viewport can be viewed in the principal-axis (PC) frame or the crystallographic (a, b, c) frame.**
  A **Frame** switch (PC | Crystal) in the plot header flips the axis triad, the shadow-box wall projections,
  the look-down camera buttons, and **Reset view** together — PC1/PC2/PC3 in PC mode, a/b/c in crystal mode. In
  crystal mode the box + walls switch to an orthonormal frame built from the unit cell and the SAME 3D KDE
  density is re-binned onto the a/b/c planes (`projectDensityOntoFrame`), so the wall shadows are the honest
  crystal-plane projections of the displayed density; the a/b/c rods and the look-down-a/b/c camera use the true
  cell edges. Reset view frames whichever box is active corner-on, and switching frames snaps to that frame's
  default view. a/b/c are keyed by rod color in the viewport legend (no 3D letter labels). The crystal-frame directions come from `src/pcaCrystalFrame.js`
  (`unitCellVectors`), which derives the unit-cell vectors in the shared Cartesian basis the PCA axes and
  density already live in (`src/__tests__/pcaCrystalFrame.test.js`).
- **The Site ellipsoids panel gains crystallographic axes + reset/save.** An opt-in a/b/c gizmo at the unit-cell
  origin (its own **a b c** toggle, keyed by color in the panel legend) and **Reset view** + **Save** (PNG,
  1× / 3×) controls mirroring the main panel.

Old coordinates-only `.rmc6f` files reconstruct thermal-ellipsoid sites; PCA controls regrouped with a shell contrast knob.

- **The oldest `.rmc6f` files (coordinates only) now drive PCA-KDE and Atomic Density.** That format drops the
  reference-site and per-atom cell columns entirely (`id element [label] x y z`, 6 fields), so there was no
  per-site grouping — the file failed to load at all. `parseAtomLine` now also accepts those short lines, and
  when a file carries no reference numbers the PCA parser reconstructs sites by folding every atom into one unit
  cell and clustering per element (periodic, full cell-metric minimum image, circular-mean unwrap of each
  cluster). Each reconstructed site reports its member count against the one-per-cell expectation from the
  supercell: `27/27` is a clean crystallographic site, while `162/27` flags atoms that do not resolve into
  separate sites at the chosen distance — close sites, or an orientationally-disordered group such as a rotor
  shell, whose "ellipsoid" is a shell best read from the KDE. A **Cluster** distance knob (shown only for such
  files) tunes the grouping; the site list and a badge surface the count and a clean / merged-or-disordered
  label, and the default selection prefers a clean site. Atomic Density needed only the coordinate fold, which
  no longer depends on the cell columns. New synthetic, oblique-cell, and real SF6 190 K (rotor-phase) cases in
  `src/__tests__/rmc6f.test.js`. (The Flask `/api/pca` backend does not reconstruct; the browser worker path does.)
- **PCA Ellipsoid controls regrouped for a clearer layout.** The wireframe and the density painted on its
  surface now share one **Ellipsoid** group (Wireframe · Level · Color · Shell · Colormap · Contrast), with
  every toggle labeled by what it toggles. A new **Contrast** knob stretches the KDE-shell colormap around its
  mid-tone — a single control over the effective vmin/vmax — so faint departures from the harmonic ellipsoid
  stand out. Turning on the **Isosurface** now also clears the wireframe and shell for a clean volume view.
- **The Displacement-statistics panel no longer clips.** It now sizes to its content (so the full covariance
  matrix, all three principal axes, and the summary line always show), and the Site-ellipsoids 3D view below it
  re-fits to the remaining space even though its PCA arrives asynchronously from a worker.

Dashboard plots labeled by the fit function declared in the run-control `.dat`.

- The RMCProfile run-control `<stem>.dat` records the correction / fit-function form per dataset
  (`> DATA_TYPE :: G(r)` with `> FIT_TYPE :: D(r)` means the data is fit as D(r)). The dashboard now uses
  that fit type as the plot's heading and y-axis label, so a `.gr` file fit as D(r) reads **D(r)**, not
  G(r) (F(Q), S(Q), … likewise). `parseRunSettings` already extracted it; a new `fitTypeByFilename` maps it
  by file name and the run assembly pairs it onto the matching plot file.
- Finding the right `.dat` among the many in a run folder (chi2.dat, optimization.dat, …) uses the existing
  stem match (`<rmc6f-stem>.dat`) first, then a capped content scan of the other `.dat` files as a fallback,
  reading only each file's head so a large data `.dat` is never read in full.
- Any `.gr` / `.sq` / `.fq` STOG data file now loads as a dashboard plot, not just the default `scale_ft.*`
  names — runs commonly use descriptive data-file names (e.g. `PMN_300k_rmc_..._v2.gr`).

Robust `.rmc6f` parsing, dashboard box-zoom, Atomic Density render fix, and PCA polish.

- **Older `.rmc6f` files now load.** The only structural difference from the current format is the per-atom
  type label (`id element [type] x y z ref cx cy cz`, 10 fields); older 2018-era files omit it (9 fields).
  Atom-line parsing now indexes the reference number, cell indices, and coordinates from the END of the line,
  tolerating any number of label columns. A new shared `web_app/frontend/src/rmc6f.js` (`parseAtomLine`) backs
  both the structure and PCA browser parsers; Python `iter_rmc6f_atoms` (used by `pca_kde` and `kde`) gets the
  same treatment. Covered by `src/__tests__/rmc6f.test.js` and `OldFormatRmc6fTests` in `tests/test_parsers.py`.
- **Dashboard charts support box zoom.** Drag a rectangle to zoom into that region on both axes (a thin
  horizontal or vertical drag still zooms just that axis); series are clipped to the plot area, and Reset
  zoom / double-click restore the full view.
- **Atomic Density first-visit render fix.** The KDE-slice / slab canvases and the 3D model measured their size
  once on mount with no observer, so a first visit before layout settled could leave the 2D panels below
  display resolution and the 3D model at 0×0. They now re-measure via a `ResizeObserver` (2D redraw; 3D updates
  the camera + renderer in place), matching the pattern the PCA viewport already uses.
- **PCA Ellipsoid polish:** a **Black** wireframe option; the wireframe is drawn opaque while the KDE shell is
  on (crisp cage over the colored surface); the **KDE shell** toggle's `?` moved out of the switch `<label>`
  so the toggle click is no longer swallowed by the help badge, and its popover opens leftward (`align="end"`);
  Site-selector markers are drawn as each site's calculated thermal ellipsoid, with the selection glow + triad
  scaled to the marker so a soft site can't outgrow them; and the U_iso / B_iso subscripts render correctly
  (the stat label was an inline-flex box that dropped `vertical-align` on `<sub>`).

PCA Ellipsoid: selectable / KDE-projected ellipsoid, Site selector; taller dashboard charts.

- **KDE shell** (new, optional): a "KDE shell" toggle paints the KDE density onto the p% ellipsoid surface
  (per-vertex trilinear sample of the volume, current colormap), so anharmonic departures from the harmonic
  ellipsoid read as hot/cold patches on the shell. The coloring is stretched to the shell's own density range
  (not the global 0..vmax), so the small variation across a near-iso-probability shell reads at full contrast
  (a Gaussian cloud stays near-flat, as it should). It shows the same density as the isosurface from the
  outside and would occlude it, so the two are mutually exclusive — enabling one switches the other off, and
  both-off is allowed (wireframe + wall projections only).
- **Ellipsoid wireframe color** is now chosen from the controls bar (swatch + dropdown: Amber, White, Silver,
  Cyan, Violet), kept transparent so it never crowds the isosurface. The default is Amber (`#ff7a1a`) — the
  same warm tone as the selected-site highlight in the Site selector, so the cage reads as "this site."
- The **Unit cell** side panel is renamed **Site selector**, since its job is picking the reference site to
  analyze (click an atom to load its PCA-KDE).
- Dashboard charts fill the panel width on wide (≥ 1500 px, 16:9) screens instead of being capped short and
  letterboxed with empty side margins: the 8:5 plot's height cap goes 360 px → 540 px, so it grows to fill the
  card at the same aspect (e.g. 1080p: 576×360 with side gaps → 618×386 filling). Fills fully through 1440p.

PCA Ellipsoid statistics detail and a 16:9 layout; Demo button beside the page tabs.

- The **Displacement statistics** panel now shows the full covariance tensor U (Å², Cartesian x/y/z with the
  variance diagonal highlighted) and a **Principal axes** table — each PCA component (eigenvector) as a row
  in Cartesian x/y/z with its eigenvalue λ (Å²), RMS amplitude (Å), and per-axis excess kurtosis κ. PC1/PC2/
  PC3 rows are color-keyed to the 3D triad and the viewport legend. The old single-line "RMS axes" and
  "Eigenvalues" rows fold into that table (now per-PC columns); U_iso / B_iso / Anisotropy / Non-Gaussianity
  stay as the scalar summary. Both the Flask and browser-worker paths already return `covariance` / `axes` /
  `eigenvalues` / `rms` / `excessKurtosis`, so the detail appears in every runtime mode.
- The page now fills the screen height instead of sprawling wide. The 3D viewport takes three-quarters of the
  width on the left; the statistics and unit-cell panels stack as two **equal-height** rows in the remaining
  quarter, with the viewport spanning both so all four outer edges line up. The grid flexes to fill the page's
  leftover height, so the view roughly fills a 16:9 screen. The four scalar figures sit in a 2×2 grid and the
  statistics body scrolls within its row if a small window can't show every table; below 980 px it stacks into
  one column at natural heights and the page scrolls.
- The header **Demo** toggle moves into the left cluster just after the AI Assistant tab (beside the page
  tabs), leaving Live Data + the run-folder field as the right-hand data-source group.

PCA Ellipsoid: reset the main view on a new dataset, and a deterministic camera reset.

- Loading a different run now returns the main 3D panel to the default body-diagonal view instead of
  inheriting the previous model's orbit. The reframe is keyed on the run's identity (App's `runId`),
  so it fires even when the new run's first site shares a reference number with the old one, while a
  Live Data refresh of the same run and any slider/layer tweak leave the view untouched.
- Root-caused and fixed why the reset sometimes landed rotated off the default: `frameMainCamera` set
  the pose and then called `OrbitControls.update()` once, which applied the control's leftover damped
  orbit velocity to the new pose. It now flushes that velocity first (a damping-off `update()` that
  zeroes the pending delta), so the reset is deterministic regardless of prior momentum — the same
  reason it also fixes the "Reset view" button landing slightly off right after a drag. This let the
  reframe live in one place (the scene rebuild) and removed the visibility-tracking scaffolding.

PCA displacements in the AI Assistant context.

- The run context sent to the local model now carries a `pca_displacements` section: one entry per
  reference site from the PCA Ellipsoid analysis (isotropic ADP `U_iso_A2`, the three principal RMS
  amplitudes `rms_axes_A`, `anisotropy`, and mean-excess-kurtosis `non_gaussianity`), ranked most
  non-Gaussian first so the anharmonic / split sites lead and survive budget trimming. This is
  information the model previously couldn't see — symmetry gives only mean displacement per Wyckoff
  orbit, not the anisotropy or the anharmonicity. The system prompt gains a matching bullet so the
  model interprets the kurtosis correctly, and the character-budget ladder trims the least-anharmonic
  PCA sites (then their explanatory note) alongside the other evidence.
- Wiring: the PCA Ellipsoid page publishes its computed per-site table upward (`onSitesChange`), App
  holds it, and the AI Assistant page threads it into `buildRunContext`. The section is present once
  the PCA Ellipsoid page has been opened for the run (its analysis runs there), and follows the active
  dataset. Covered by four new cases in `src/llm/__tests__/runContext.test.js`.

PCA Ellipsoid page refinements.

- The main 3D view has a **Reset view** button in the panel header that returns the camera to the
  default body-diagonal framing after you orbit or zoom away — the same framing applied automatically
  when the site changes, now available on demand.
- Fixed the PCA Ellipsoid not recomputing for a newly loaded dataset. The static-mode worker cached
  parsed displacement clouds by file *path*, so a different run reused the previous model's clouds. The
  cache is now content-addressed (keyed on a signature of the `.rmc6f` text), and the page ties the
  loaded text to the file it came from, so requests never run against the previous dataset. Covered by
  `src/workers/__tests__/pcaKdeWorker.test.js`.

PCA / thermal-ellipsoid KDE computation engine.

- New engine turns per-site RMC displacement clouds into anisotropic displacement tensors (the
  thermal ellipsoids) and a smooth 3D probability density. An atom's offset from the average
  structure is `coords − cellIndices/supercell` folded over the supercell boundary, mean-subtracted
  per reference site, then mapped to Cartesian Å; the covariance of each site's cloud is its ADP,
  and a Gaussian KDE of the cloud is the density the ellipsoid only approximates. PCA and KDE are
  standard tools; the specific per-site analysis and the shadow-box visualization below follow
  **Maksim Eremenko's PCA_KDE utilities**
  (<https://github.com/MaximEremenko/Utilities/tree/main/RMCProfileUtilities/PCA_KDE>) — we followed
  that approach and reimplemented it independently into the toolkit's dual-mode architecture, not a
  port of his code (his `KDE.js` evaluates a full multivariate Gaussian KDE; ours factorizes it).
- **The KDE is separable, and exact.** Sampling on a grid aligned with the cloud's principal axes
  makes SciPy's `H = factor²·C` bandwidth diagonal, so the 3D Gaussian factorizes into three 1D
  kernels and the volume becomes their tensor product: `N·3·grid` exponentials contracted through
  BLAS instead of `N·grid³` kernel evaluations. Measured **36–51× faster than a naïve
  `scipy.stats.gaussian_kde` sweep of the same volume, with max abs error ~1e-13** (machine
  precision) — it is the same estimator, not an approximation. A 52-site, 52 000-atom configuration
  builds all clouds in ~80 ms and all 52 ellipsoids in ~1 ms; a single-site 48³ volume solves in
  ~7 ms server-side.
- `rmc_toolkits/pca_kde.py` is the source of truth: `load_site_displacements`, `site_ellipsoids`
  (batched, one pass), `pca_kde_volume` / `site_pca_kde` (volume + PC-plane projections + iso
  thresholds by enclosed probability mass and by raw density). Axes are sign- and
  handedness-canonicalized for reproducibility; a floored eigenvalue keeps flat/linear clouds
  non-singular and flags them `degenerate`.
- `web_app/frontend/src/workers/pcaKde.js` is a straight JS port for static mode (a 3×3 Jacobi
  eigensolver stands in for `eigh`), driven by `pcaKdeWorker.js` off the main thread — ~96 ms for a
  1000-point 48³ volume in-browser, no GPU needed. `tests/test_pca_kde.py` and
  `src/workers/__tests__/pcaKde.test.js` both assert the separable volume equals the full
  brute-force estimator to round-off (15 + 14 tests), plus site extraction, supercell-boundary
  folding, mass normalization, and known-anisotropy recovery.
- Flask exposes `/api/pca/sites` (per-site ellipsoid table) and `/api/pca/kde` (one site's or one
  element's volume), both behind a per-(path, mtime) LRU cache.
- **PCA Ellipsoid page** (new top-level tab; the "KDE / 3D" tab is renamed "Atomic Density"). A site picker over all reference sites, an
  ellipsoid summary table, a Three.js scene showing the KDE isosurface (extracted by a self-written
  marching-cubes module — `workers/marchingCubes.js`, since Three's `MarchingCubes` only builds
  metaballs) nested with the p% thermal-ellipsoid wireframe and PC-axis triad, and the three
  PC-plane KDE projections as heatmaps. An isosurface-mass slider sweeps the enclosed-probability
  threshold. Works in both runtimes: Flask endpoints in backend mode, the `pcaKdeWorker` +
  `pcaKde.js` engine off the main thread in static mode (verified on the bundled Demo run).
- **Non-Gaussianity (excess kurtosis)** is reported per site. It quantifies why a KDE isosurface can
  sit inside its harmonic ellipsoid: the covariance is inflated by fat tails while the isosurface
  tracks the peaked core. Verified on `data/5K_try1` — Nb sites read ~5–7 (strongly anharmonic,
  isosurface well inside the ellipsoid), near-isotropic Ga sites read ~1 (isosurface ≈ ellipsoid),
  and a synthetic Gaussian cloud reads ~0. Marching cubes has its own sphere/normals tests; the
  isosurface-vs-ellipsoid scaling was validated to agree with theory (KDE surface ≈ 1.06× the
  ellipsoid for a true Gaussian).

Periodic boundary conditions for the KDE slice.

- The KDE under-counted near cell boundaries: folding atoms into one unit cell drops their
  periodic neighbors, so the density decayed toward faces/edges/corners (an edge blob showed
  roughly half the interior amplitude, a corner roughly a quarter), and a slab centered on the
  z=0 face missed the atoms just below z=1 entirely.
- Both KDE paths now tile periodic images from the 26 neighbor cells within a margin that covers
  the kernel reach and the slab depth (`min(0.5, max(0.1, 2*bw, thickness))`). In
  `rmc_toolkits/kde.py`, `_augment_periodic_images` feeds `oriented_kde_slice`, and `kde_slice`
  reports `slabCount` as unique source atoms and rescales the density (SciPy's `gaussian_kde`
  divides by every fit point, images included). In `localKdeWorker.js` the same augmentation runs
  before `makeSlab`, and the image factor rides in the kernel normalizer, so the WebGPU path
  inherits the fix unchanged. `kde_slice` called directly (no `source_index`) keeps the old
  truncated behavior for generic point clouds.
- Verified on `snao.rmc6f` at a z=0 slab: slab population 3110 → 6034 atoms, corner density
  20.8 → 55.4, edge-column peaks now match interior peaks (~73 vs ~75), opposite edges agree to
  ~5% (fit-subsampling noise). New tests: in-plane wrap symmetry and depth-wrap slab selection in
  `tests/test_kde.py`, plus a mirrored vitest suite for the worker (whose `onmessage` registration
  is now guarded so tests can import the module outside a worker).

Source data file names on figures (user-experience feedback).

- Every dashboard plot card shows its source file name in monospace under the title — the two
  otherwise-identical `EXAFS Q-space` cards are now distinguishable — with the full path on hover.
- The collapsed R-value strip shows its source log name(s); a combined multi-log strip lists every
  parsed log.
- The Model information title cell now displays the source `.rmc6f` file name (previously hidden
  in a tooltip), which also covers the KDE / 3D panels since they all derive from that one file.
- File names sit outside the card `h3`, so "Save all figures" export names are unchanged.

Math rendering in AI Assistant replies.

- Assistant Markdown now typesets LaTeX: integrated `remark-math` + `rehype-katex` so
  dollar-delimited math (e.g. `$R_{wp}$`, `$\chi^2$`, `$\text{\AA}$`) in replies renders as real
  math instead of literal text — including inline math inside GFM table cells, which is how models
  tend to format result summaries. KaTeX runs *after* `rehype-sanitize` (raw → sanitize → katex) so
  model HTML stays fully sanitized while KaTeX emits its own trusted markup, and KaTeX fonts are
  bundled locally so the static GitHub Pages build needs no CDN.

Demo run and refreshed screenshots.

- Added a header **Demo** toggle that loads a bundled GaTa4Se8 250 K example run (under
  `web_app/frontend/public/demo/`) so first-time visitors see a populated dashboard; a second click
  clears it.
- Recaptured the Dashboard and KDE/3D screenshots against the demo run and added an AI Assistant demo
  GIF (a local Ollama model summarizing the run as a table with inline LaTeX math). Renamed the screenshots to
  `assets/rmc-toolkits-dashboard-demo.png` / `assets/rmc-toolkits-kde-demo.png` (new paths bust the
  stale README image cache) and removed the unused legacy screenshots from `assets/`.

Relicensed to AGPLv3.

- Switched the project license from MIT to the **GNU Affero General Public License v3.0** to keep
  the code open while requiring that modified versions — including those run as a hosted network
  service (AGPL §13) — release their source. Updated `LICENSE`, `pyproject.toml` (license field +
  trove classifier), the README badge/notice, and the in-app footers.
- Added an **About & documentation** link in the Dashboard and KDE/3D footers pointing to the repo,
  which also surfaces the source for AGPL §13 network-use compliance.

RMCProfile-first positioning.

- README, QuickStart, roadmap, and agent notes now frame the current app as an RMCProfile
  modeling-output dashboard. STOG is treated as deferred preprocessing/legacy parser support, not
  a current RMC modeling feature.
- The dashboard excludes STOG plot kinds from the visible plot list and assistant context while
  keeping the underlying parser support in place.

RMCProfile wording cleanup.

- Replaced RMCProfile "refinement" language in docs and assistant prompts with RMC modeling /
  atomistic configuration optimization wording, to avoid implying Rietveld-style parameter
  refinement.
- Renamed the assistant's move-counter context from `refinement` to
  `configuration_optimization`.

Gentler AI Assistant startup.

- The automatic on-load connection probe no longer shows a red "Connection failed" before the user
  has set anything up. The status reads a neutral **"Connect a model to start"** (grey dot), and the
  settings drawer stays clean; the full error + setup hint appear only after the user actively
  presses **Test**.

Assistant "Beta" badge and a cleaner chat speaker label.

- The header badge now reads **Beta** (was "Experimental") in **amber** — the complement of the
  app's cobalt accent, and the same palette as the cloud-provider notices — so the beta reminder
  stands out instead of blending into the blue theme.
- Chat messages are labeled **Assistant** rather than the raw model id (which stays visible in the
  header's model switcher), for a cleaner transcript.

Run history and run-control settings in the AI Assistant's context.

- **`configuration_optimization` block** from the `.rmc6f` header: moves generated/tried/accepted
  plus the derived **acceptance ratio** and **accepted moves per atom** — the standard gauge of
  whether the configuration has been sampled long enough — and accumulated running time.
- **`run_settings` block** from the RMCProfile run-control file. The correct `.dat` is the one whose
  basename matches the chosen structure stem (so `chi2.dat`, `optimization.dat`, … are never picked
  up). Extracts title/material/phase/temperature, **minimum distances labeled per element pair**
  (hard closest-approach constraints — the model is told a g(r) peak pinned there may be
  constraint-limited), max move sizes per element, time/save limits, flags, and the fitted-data
  blocks. Static mode only, like the rest of the browser-parsed context.
- System prompt teaches both blocks; 9 new tests (header counters, `.dat` parser, stem-matched
  selection, derived stats, pair labeling).

Chat rendering fixes and a window-filling chat box.

- Markdown now parses **and sanitizes** inline HTML (`rehype-raw` + `rehype-sanitize`), so `<br>`
  line breaks inside table cells render instead of collapsing onto one line; scripts, event
  handlers, and images are still stripped (no injection surface, no external loads from model
  output) — covered by tests.
- The AI Assistant chat box now **fills the window height** (message log scrolls internally, the
  composer stays pinned at the bottom) and the column is a bit wider (880 → 1000 px).

Rendered Markdown in AI Assistant chat replies.

- Assistant answers now render **Markdown** — GitHub-flavored **tables**, lists, code blocks, inline
  code, headings — via `react-markdown` + `remark-gfm`, so a "summarize this as a table" request
  shows a real table instead of raw `| --- |` text. Rendering builds a React tree with raw HTML
  disallowed, so the no-injection property is preserved (covered by a test); user messages stay
  plain. (eslint: enable `ignoreRestSiblings` for the `{ node, ...props }` component overrides.)

App identity — a minimalist pair-distribution "wave" mark.

- Replaced the placeholder "R" brand with a single shared g(r)-wave SVG (flat brand blue, `#2563eb`):
  a base-aware favicon (`public/favicon.svg`), the header brand mark, and the AI Assistant
  avatar / empty-state icon. Removed the default `vite.svg`.

ChatGPT-style "Thinking" indicator for reasoning models.

- Reasoning models (e.g. qwen3 via Ollama) stream their chain-of-thought in a separate field; the
  streaming client now yields structured content/reasoning chunks and the chat shows a collapsible,
  shimmering **"Thinking"** panel (live reasoning) that collapses to **"Thought for Ns"**
  (re-expandable), plus a bouncing-dots indicator for non-reasoning models slow to their first token.
  Reasoning is preserved on each turn.

Local-distortion evidence in the AI Assistant's run context.

- The context JSON gains a **`symmetry` block** — space group, the symmetry-vs-tolerance **ladder**
  (distortion magnitude + character), and the Wyckoff-orbit **sites ranked by rms displacement** —
  plus **`pair_correlations`**: each partial-PDF pair's first g(r) peaks next to the
  nearest-neighbour distance the average structure predicts, so the model can point at the sites
  and pairs participating in short-range correlations.
- Per-site rms displacements (`dispA`) fall out of the circular-mean accumulators
  `structureFromRmc6f` already had — derived per axis from the resultant length, scaled to Å; the
  symmetry orbits now carry `members` (basis indices) so the context can aggregate per orbit.
- New `llm/context/pairCorrelations.js` (average-structure NN distances over periodic images +
  a smoothed local-maxima peak finder); context budget raised to ~4.5k chars with an explicit trim
  order (extra peaks → ladder middle rungs → history → low-displacement sites → datasets, each
  recorded with an `*_omitted` count); system prompt teaches the model how to read the new blocks.
- Tests: 51 passing (peak finder, NN distances, displacement math on synthetic .rmc6f fixtures,
  symmetry-block assembly, trim order).

AI Assistant rework (shipped as PRs #3/#4, catching the log up): moved from a dashboard card to a
dedicated **AI Assistant page** (chat-only — Summary/Report tabs removed, everything happens in the
chat), modern chat UI with a persistent connection bar (status dot, model switcher, settings gear,
auto-connect), redesigned Connection Settings, and optional **cloud providers** (OpenAI, Gemini)
with Bearer API-key auth and explicit data-leaves-your-device warnings; local Ollama/LM Studio
remains the default. Ollama CORS + Safari setup documented in-app and in QuickStart.

Experimental AI Assistant — local LLM in the data-monitor pipeline.

- Added `web_app/frontend/src/llm/`, a self-contained experimental module connecting the dashboard
  to a **local LLM** (Ollama `/v1` or LM Studio, both OpenAI-compatible) **directly from the
  browser** — no server, no API keys; run-derived summaries go only to the model server the user
  runs. Four features: one-click **run summary/assessment** (streamed), **chat Q&A** with the run
  context injected, one-click **Markdown run report** (deterministic metrics tables + clearly
  labeled AI narrative), and a **live convergence watchdog** badge on the R-value card
  (slope heuristics are the source of truth; the LLM only writes the note; piggybacks on the
  existing Live Data poll).
- Pipeline seams are explicit as a learning artifact (see the module README): context builder with
  history downsampling + a ~3k-char budget → prompt templates → hand-rolled SSE streaming client →
  streaming UI with a "context sent to the model" inspector. Connection test translates failures
  into actionable CORS/setup hints.
- Dashboard mounts the collapsible **AI Assistant** card (zero network activity until used) and the
  watchdog badge; strict import boundary (props in, one `figureExport` helper out) keeps the module
  extractable to its own repository.
- Added **vitest** (first frontend unit tests, 43 passing) covering the context builder,
  convergence heuristics, SSE parsing, prompt templates, and report assembly; CI now runs
  `npm test` between lint and build.

## v0.3.0 — 2026-07-02

Symmetry analysis in the browser.

- Added a **Detected SG** card beside *Model information* on both the Dashboard and KDE/3D pages: a
  client-side, table-free **FINDSYM-like space-group finder** (ported from the RMC-phonon-dynamics
  project — no spglib/WASM). From the folded `.rmc6f` basis it detects the space group (H–M symbol +
  number), point group, and operation count, and draws an interactive **tolerance ladder** — each
  brick is the space group that holds over a range of atomic-position tolerance; click a rung to
  select it. New frontend modules `symmetry.js` (the finder) and `symmetryModel.js` (structure →
  conventional cell + basis glue); `browserData.structureFromRmc6f` now also returns a
  circular-mean per-reference-number basis.
- The space-group **number table is complete for every producible symbol** (all point-group ×
  allowed-centering combinations), and the headline always shows the point group even when a number
  is unavailable.
- The **Detected SG tolerance selection persists** when switching between the Dashboard and KDE/3D
  pages (shared via context).
- KDE/3D: **contours and log-scale density are on by default.**

## v0.2.0 — 2026-06-26

Figure export, axis labels, and browser-first positioning.

- Added per-figure **Save** controls. Dashboard charts (inline SVG) export as **PNG** (raster) or
  true-vector **SVG**; the KDE/3D panels (canvas + WebGL) export **PNG** at native or **3×**
  resolution. The earlier raster-PDF attempt was dropped — a PDF of a chart isn't true vector, and
  SVG is the proper vector format (convertible to PDF offline if needed).
- Added **Save all figures** in the dashboard's *Loaded N plot files* header. It bundles every
  visible chart into a single `.zip` (PNG or SVG), avoiding the browser's multi-download blocking.
- Added supporting modules: `figureExport.js` (SVG/canvas rasterization + standalone-SVG
  serialization with inlined computed styles), a shared `SaveMenu` component, and a dependency-free
  store-method ZIP writer (`zipArchive.js`, with CRC-32).
- High-resolution panel export re-renders the same drawing onto a 3× offscreen canvas (KDE slice,
  slab) or re-renders the Three.js scene at a higher pixel ratio (`preserveDrawingBuffer`), so the
  output is genuinely higher-resolution rather than upscaled.
- Gave proper vertical-axis labels to plots that previously showed a generic `data`: PDF/partials →
  `G(r)`, x-ray/neutron S(Q) → `S(Q)`, Bragg → `Intensity`. Applied in both the Flask
  `/api/plot/data` endpoint and the browser static-mode parser.
- Reframed the README and roadmap around the hosted browser app (GitHub Pages) as the primary,
  no-install way to use the dashboard; the Flask backend and Python package are now positioned as
  optional local/advanced paths. Added a Phase 7 plan for interactive STOG reduction and guided
  RMCProfile run setup.
- Capitalized **RMCProfile** consistently in the app title, header, and messages.

## v0.1.0 — 2026-06-21 (first tagged release)

First public, releasable version: the reusable Python package, Flask API, interactive dashboard,
server-side + browser (WebGPU/CPU) KDE, and Three.js structure viewer with a draggable slice band.

RMCProfile EXAFS dataset output parsing and large-structure display (added late in the v0.1.0 cycle):

- Added plot-kind detection for RMCProfile EXAFS dataset output files (`*-EXAFS-*_Q_OUTPUT.csv` and
  `*-EXAFS-*_R_OUTPUT.csv`) in both the Python package and browser static mode.
- Added `read_exafs_csv`, which handles Q-output files with a descriptive title row before the column
  header and R-output files with real/imaginary/modulus transform columns.
- Added backend and frontend chart labels for these dataset outputs:
  - Q-space: `k (Å^{-1})` vs `χ(k) k²`.
  - R-space: `r (Å)` vs `FT[χ(k) k²]`.
- Added dashboard file badges/order and tests for parser, plot metadata, file listing, and
  `/api/plot/data`.
- Raised the structure endpoint and Structure page point limits to 1,000,000 and raised the Slab In
  Cell canvas draw cap to match, so moderate structures such as `data/RMC/snao.rmc6f` render all
  returned atoms instead of every other point.

Release engineering added in this version:

- **License:** added an MIT `LICENSE`.
- **Packaging:** added `pyproject.toml` so `rmc_toolkits` is pip-installable (`pip install -e .`) and
  exposes `rmc_toolkits.__version__`.
- **CI:** added `.github/workflows/tests.yml` — runs the Python test suite plus frontend lint/build
  on every push and PR to `main`.
- **Tests:** sample-data-backed tests now skip cleanly when the gitignored GNSe example dataset is
  absent, so the suite is green on a fresh clone and in CI (17 run, 16 skip without the sample).

## 2026-06-21 — Draggable slice band + structure-view color work

- **Slab In Cell drag-to-move slice.** Users can now grab the highlighted band in the *Slab In Cell*
  panel and drag it along the slice axis to set the slice position (`zCenter`) live; thickness is
  unchanged. Implemented by adding an `invert()` to `makePlaneMapper` (recovers plane coordinates
  from a screen point), publishing the band geometry each render into `slabGeometryRef`, and wiring
  pointer handlers on the slab canvas (`grab`/`grabbing` cursors, pointer capture,
  `touch-action: none`).
- **Bug fix:** the drag listeners initially never attached — the effect used `[]` deps but the slab
  canvas only renders after a structure loads (`{structure && ...}`), so `slabCanvasRef` was null at
  mount. Fixed by depending on `[structure]`. Verified end-to-end against `data/5K_try1`: drag moved
  SLICE 0.39 → 0.709 with the band and in-slab atoms updating.
- Earlier in the day: gave each element a distinct color across the slab and 3D views, added a shared
  atom color legend above the structure views, and relocated the atom badges.

## 2026-06-19 — WebGPU-accelerated browser KDE

- Moved the static-mode density-map hot loop (`O(grid² · samples)`, up to ~400M `exp()` per slice)
  onto the GPU. Each grid cell is independent → one WebGPU compute-shader invocation per cell.
- Added `workers/gpuKde.js`: inline WGSL shader, lazily-cached adapter/device init, `computeDensityGpu`
  (writes buffers, dispatches, reads back via `mapAsync`, reshapes to the same `density[y][x]` grid as
  the CPU path), and a `shouldUseGpu` work-size heuristic.
- Refactored `localKdeWorker.js` into `computeDensityCpu` + a GPU-or-CPU branch; made `computeKde`
  and `onmessage` async. Same `{ id, result }` message shape, so `StructurePage.jsx` needed no
  changes (its request-`id` guard already drops stale replies). Result carries `backend: 'gpu'|'cpu'`.
- GPU used only when `grid·grid·samples >= 2_000_000`; init attempted once per worker and cached.
- **Robustness is the design point:** missing `navigator.gpu`, no adapter, device/shader error, lost
  device, or sub-threshold work all fall back to the JS loop with identical output; no GPU failure
  surfaces through `worker.onerror`.
- Verified (Chromium/Metal): GPU(f32) vs CPU(f64) parity ~1.7e-6 relative across grids 16/120/260;
  density step ~59× faster at grid 120 and ~107× faster at grid 260. Still subsamples to 6000 slab
  points and uses f32 — a visualization path, not a substitute for server-side SciPy KDE.

## 2026-06-18 — Live Data and loaded-file controls

- Added optional Live Data monitoring for the local Flask dashboard: polls the selected folder for
  supported-file changes and refreshes metadata/charts without a manual reload.
- Added file mtime + size to `/api/files` so the frontend can detect changes cheaply.
- Added a collapsed `Loaded N plot files` panel; expanding shows detected files as badges that can
  hide/show their chart.
- Aligned local folder-selector copy with the hosted dashboard (`Run folder` / `Select Folder`) and
  hardened the native picker startup so a missing default folder falls back to the nearest existing
  directory.

## 2026-06-16 — Hosted static dashboard

- Added the GitHub Pages workflow (`.github/workflows/pages.yml`); Pages source set to **GitHub
  Actions** so the built React/Vite app is served instead of a Jekyll README page.
- Added static-mode local file loading (`browserData.js`): open the hosted dashboard, select a local
  run folder, and parse RMCProfile CSV/log/STOG outputs with no Python and no upload.
- Extended static mode to parse uploaded `.rmc6f`, populate the model summary, and render the slab
  projection + Three.js 3D view.
- Added `workers/localKdeWorker.js`: off-thread browser KDE capped at 6000 slab points, deterministic
  pseudo-random sampling, contour segments for the existing overlay.
- Kept Flask `/api/kde/slice` as the reference server-side SciPy path.

## 2026-06-15 — Palette and docs

- Three.js atom palette → Nature-style (Ga deep blue, Nb vermillion, Se teal).
- Refreshed `assets/rmc-toolkits-KDE.png`; updated README + frontend docs for the local preview flow.

## 2026-05-27 — UI refresh and interactive dashboard

- Bright/dark theme system via CSS variables with a persisted header toggle.
- Removed the sidebar-first workflow; header now hosts data-path input, Dashboard/KDE nav, and theme
  toggle. App defaults to the `data` path.
- Added `GET /api/plot/data` (parsed series + normalized scientific labels: `χ`, `Å`, `Q (Å⁻¹)`).
- Replaced PNG cards with `InteractivePlot.jsx` (SVG, hover readouts, legend toggles, integer x
  ticks, drag-to-zoom + reset). Simplified the Dashboard to a three-card grid.
- Reworked KDE/3D into three side-by-side panels (XY slice, slab x-z projection, Three.js model);
  2D panels preserve lattice aspect ratios. Added gray slab-edge outlines in the 3D model.

## 2026-05-27 — Real server-side KDE

- Added `rmc_toolkits/kde.py`: loads unit-cell-folded cartesian (Å) positions from a `.rmc6f`
  (optional element filter) and computes an XY `scipy.stats.gaussian_kde` density for a z-slab
  (ported from `src/RMC_KDE.py`). Returns plain arrays (grid, extent, contour polylines via
  `contourpy`, slab count, vmin/vmax); fit subsampled to 6000 points for slider responsiveness.
- Added `GET /api/kde/slice`. `z`/`dz` are cell-edge fractions converted to Å internally; loaded
  positions cached per (path, mtime, element) with `lru_cache`.
- Made the backend port configurable via `RMC_TOOLKITS_PORT` (default 5000).
- Replaced the browser box-blur "density" with the fetched real KDE grid (colormap LUT + contour
  overlay); fetches debounced with `AbortController`. Added bandwidth/colormap/grid/contour/log-scale
  controls and `colormaps.js`. Default z-slice auto-snaps to the densest band on load.
- Added `.venv/`, `__pycache__/`, `*.pyc` to `.gitignore`.

## 2026-05-27 — Package + tests foundation

- Added `rmc_toolkits/parsers.py` (RMC CSV/log/STOG, `.rmc6f`, `Frac*.txt`) and
  `write_frac_from_rmc6f` conversion.
- Added `rmc_toolkits/plots.py` (plot detection, figures, metrics, PNG serialization); backend now
  calls into the package instead of duplicating logic.
- Added backend endpoints: `/api/health`, `/api/files`, `/api/plot`, `/api/plot/metadata`,
  `/api/convert/frac`, `/api/structure`, `/api/kde/slice`; plus a data-root guard
  (`RMC_TOOLKITS_DATA_ROOT`).
- Frontend uses `VITE_API_BASE_URL` (default `http://localhost:5000`); file-path field no longer
  fires a request per keystroke.
- Added the Dashboard all-plots view and the KDE/3D page (element filter, z-slice controls, density
  KDE canvas, Three.js folded unit cell with OrbitControls + translucent slab overlay, and the
  slab-in-cell x-z projection).
- Fixed structure sampling so the sample renders all 52,000 atoms and preserves all 52 reference
  sites (sampling grouped by reference number, not raw stride).
- Added a `unittest` suite under `tests/` covering CSV/log parsing, Rwp, `.rmc6f` metadata +
  conversion, structure loading, plot detection/metadata/PNG, and KDE loading/slice.
