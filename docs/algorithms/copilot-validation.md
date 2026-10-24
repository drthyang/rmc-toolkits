# AI Copilot validation: does the Copilot help?

The Copilot's own tests replay scripted answers. This protocol asks a different question: given a
configuration whose local structure is **known**, does the Copilot reach the right conclusion more
often than the Workbench without it, and does it avoid claiming order that is not there? It was
fixed in code (`web_app/frontend/validation/`) before any model run. The tool and skill review
that went with it is [copilot-review.md](copilot-review.md).

## Arms

| Arm | What answers | Code |
| --- | --- | --- |
| **A** | The default Workbench, no agent: the Detected SG card (0.2 Å), the PCA Ellipsoid table, the Displacement Directions page at its own defaults (the engine's resolution, its readouts from the raw histogram, since protocol 2026-10-19.17; [below](#the-displacement-directions-page-at-its-defaults)) and the partial g(r), read by fixed rules (`ARM_A_RULES`). The cluster and superstructure engines have no page, so A cannot use them | `armA.js` `runArmA` |
| **A+** (control) | Every Copilot engine run exhaustively with its own thresholds, joined by fixed rules; no model. The cluster and modulation engines return their own verdicts, which A+ maps to claims (since 2026-10-19), a modulation claim with the engine's polarisation as its direction (since protocol 2026-10-19.12). Where the app has a rule, A+ calls the app's function (`SWEEP_RULES`, since 2026-10-19: the split site, the cluster elements, the lower group's tolerance), and a Jahn–Teller claim's long − short is the shell catalogue's split of a joined shell (since protocol 2026-10-19.15), and a split site's direction is the direction map's peak only where the peak cell resolves it, else PC1 (`directionResolved`, 15°; since protocol 2026-10-19.17), and only where the tool calls the peak a finding (`directionTier`: ≥ 3σ, ≥ 1.2× isotropic, held at the neighbouring resolution; since protocol 2026-10-24.2); a diffuse verdict in antiphase is order between cells, a modulation claim with no q (`SWEEP_RULES.diffuseOrder`; in phase is correlated motion, never claimed; since protocol 2026-10-24.2). **B − A+ is what the language model adds**; B − A mixes that with the new engines | `armA.js` `runArmSweep` |
| **B** | The Copilot as the app assembles it: tools, skills, guardrails (refusals, the Workbench's window questions, answer checks) | `harness.js` |
| **C** | B without the method ([below](#arm-c-the-copilot-without-the-method)): no skills (no `read_skill`, no skills index), no method rules in the prompt, no method advice in a tool's description or check; the tools, refusals, the Workbench's questions, the evidence gate and the answer checks still run. The ablation | `harness.js` |
| **D** | The planted truth, scored as a verdict; it must score 1 (a test checks this) | `verdict.js` `truthVerdict` |

Arm A's rules are a competent reading, not a straw man. A site is called split only when its
ellipsoid is prolate (PC1/PC2 ≥ 1.3) **and** its PC1 kurtosis is clearly negative (≤ −0.5), so an
anisotropic Gaussian is thermal motion. A split g(r) shell counts only when the average structure
has no distances that far apart, and an off-centred centre explains its own split shell.

The live arms leave out the view tools, which need a screen. The simulated user skips every
question, so the model uses its stated default. The model is told nothing about the case: every
case is loaded as the same run, `rmc_run`, with files `rmc_run.rmc6f` and
`rmc_run_PDFpartials.csv` (`runFolder.js` `RUN_NAME`), because the run name reaches the model in
its context (`"run": "rmc_run"`) and a case id names the case's kind. The record keeps the case
id; the model never sees the record. So no truth reaches the model except through the tools
(checked by `caseBlind.test.js`, below).

**The run is the app's.** `runFolder.js` `loadRun` writes a case's files as the `rmc_run` folder
and opens it through the run module's in-process adapter (`src/openRunFolder.js`), the same module
that hands the dock its run in the app (`src/copilotRun.js`, through `useCopilotRun.js`): the
Dashboard's plot parse and order, the plot kinds from `detectPlotKind`, the structure parse, the
symmetry at the view's tolerance with its ladder, the engine tools, and an engine that runs the
Copilot worker's handler in-thread. The view is `HARNESS_VIEW`, the app's defaults: the Detected
SG card at 0.2 Å with its ladder, nothing hidden, no PCA table (the PCA Ellipsoid page not
visited), Live Data off. Arms A and A+ read the same run, with the ladder off.

**The question is the app's.** Each live run is one question through the turn module
(`src/llm/copilot/turn.js`, `runTurn`), the same module the dock runs each question through: it
builds the context, the prompt, `ask_user` and `read_skill`, the history, the Workbench's
questions and the verification (`src/llm/copilot/verification/`: the evidence gate, the answer
checks, the references and the computed verdict) from the run props. The harness (`harness.js`)
is its second adapter and passes only what is its own: the case's run, its compute-only tools
(the tools over the run props with no cards, and the run's engine tools), the simulated user, the
arm's skills (arm B) and method (arm C), the verdict's words for the gate (`verdict.js`
`VERDICT_PHRASINGS`: a kind in the VERDICT line, or an empty feature list, needs its analysis as
the same claim in prose; the app's own claim catalogue reads prose only) and a model port. `run.js` builds that port as the app
does, with `providerModel` (`src/llm/provider/model.js`: the context guard, the warmer retry), and
meters it for an API model (`spend.js` `meteredModel`).

### The Displacement Directions page at its defaults

Decided 2026-10-19 (plan item 2.5; protocol 2026-10-19.17). Arm A reads a split site's direction
map on the Displacement Directions page as the page shows it at its defaults, and takes those
defaults from the page itself: `ARM_A_RULES.DIRECTION_SETTINGS` *is* the page's
`DIRECTION_DEFAULTS` (`src/directionSettings.js`), and the request is the page's own readout
request (`directionReadoutRequest`), so no number is copied into the protocol. Since this decision
the page opens at Auto, the engine's recommended resolution (`recommendedFrequency`), and takes
every readout (the peak, its direction and σ, the map test, the ± asymmetry, the Bingham
anisotropy) from the raw histogram, at smoothing 0; the smoothing redraws the map only. That is
exactly what `analyze_displacement_directions` computes when it is given no frequency, so arm A,
A+ (which reads a split site's map through the same reading, `directionsOf`) and the Copilot read
one tiling and one set of numbers; a known-answer test checks it site by site on a dev case
(`validation/__tests__/directionsDefault.test.js`). Before it, the page opened at ν = 10 and read
its numbers off the 2× smoothed map, whose peak cell, direction and σ follow the smoothing (the
engines pick the peak cell on the smoothed map).

What it changes in a reading: the reading records the resolution it read (`frequency`) and the
peak cell's radius (`radius_deg`, the "± R°" the page prints beside the direction), and a peak
that reaches 3σ on the engine's coarser tiling names the claim's direction where PC1 named it
before. A tiling names a direction only to within its cell (the engines' `peakCellRadiusDeg`, the
cell's circumscribed radius): below 504 copies Auto is ν = 1, twelve cells 37.4° in radius whose
centres lie 31.7° from every cube axis; at ν = 2 a cell is 19.3–20.9°, still wider than the 15°
direction tolerance. Arm A takes a significant peak whatever its radius, as a user reading the
page. A+ applies the tool's bound (`SWEEP_RULES.directionResolved`: `LOW_INDEX_CONE_DEG`, 15°)
and takes PC1 where the peak cell is wider, since a cell cannot localise a direction more finely
than its own size. Since protocol 2026-10-24.2 A+ also follows the tool's tiers
(`SWEEP_RULES.directionTier`, the displacement-anisotropy interview): a candidate (≥ 3σ, ≥ 1.2×
isotropic) is read again at the neighbouring resolution (ν − 1, ν + 1 at ν = 1), and the peak names
the claim only as a finding, held there at ≥ 3σ on the same axis; a hint or nothing takes PC1. The
reading records the tier, the neighbouring resolution and its σ. On dev it moves no verdict or
score: every split site's peak is a finding, and every one's cell is wider than 15° (PC1 already).
Arm A reads the page as before. The page's map opens unsmoothed, which no reading or score depends on.
The moved dev records are re-recorded in `arm-a-dev.json` (the others as recorded) and pinned in
`arms.test.js`; the tables under [Arm A results](#arm-a-results-no-model-2026-10-14-commit-b4fe09f)
were recorded under the earlier default, and the test set's records move at its next scoring.

### Arm C: the Copilot without the method

Decided 2026-10-19 (plan item 1.8; protocol 2026-10-19.10). The method is what the skills teach,
wherever the model would read it: which analysis answers which question, where an analysis takes a
parameter from (a window from `find_gr_peaks`), which window is the right one, how to judge a
result. Before this, arm C dropped only the skills, over a prompt and tool descriptions that
carried most of their rules, so it measured the skills' text, not the method. Arm C now removes the
whole method, structurally: one switch in the turn module (`config.method: 'off'`,
`src/llm/copilot/method.js`; `planTurn` decides every part from it), never text edits in the
harness. `ARMS.C` sets it; arm B and the dock follow `'full'`, byte for byte as before (their
snapshot files pass unchanged).
- **The skills:** no `read_skill` and no skills' index (as before).
- **The prompt:** no `EVIDENCE_RULES` (which analysis each conclusion needs). `TOOL_RULES` loses its
  three method sentences: the significances are calibrated one-sided σ with ≥ 3 significant; a g(r)
  peak cannot tell a distortion, so run `analyze_clusters` or `compute_bond_angles` (and
  `compare_fit` for where the fit fails); a z-score with a tiny change is not a finding. `ASK_RULES`
  loses its one ("Settle what the data can first: e.g. run find_gr_peaks, then offer its shells as
  the options."). The sentences on working with the tools (the four steps, reading `checks`, never
  repeating a call), asking, quoting by reference and the Outcome line stay.
- **The domain rules** (`SYSTEM_PROMPT`'s ground rules, shared with the watchdog; decided by the
  manager on the user's delegation, 2026-10-19): every part of every sentence is marked in
  `src/llm/prompts/system.js` (`GROUND_RULES`) as a **definition** (what a context field is: the
  quantity, its units, how it is computed, its limits as a measurement), an **inference** (which
  field to use for a conclusion, or what a pattern is evidence of) or **conduct** (how to answer:
  numbers from the context or a tool, no invented uncertainty, say what is missing, be cautious). A
  sentence that does both is split into its parts. Arm C reads the definitions and the conduct;
  every inference part goes (an inference part that ends a sentence leaves its full stop). Of the 26
  ground-rule sentences, 7 are definitions, 2 inference, 4 conduct and 13 split; for example "first
  shell gives the shell's mean distance (from r^2 g, unbiased) …" stays, and "Compare
  first_shell.mean_d_A, not the peak, with the average distance" and "A split shell … is evidence of
  short-range correlations …" go. `systemPrompt.test.js` pins the classification sentence by
  sentence. The tool-free prompt follows the method too (a model that refuses tools in arm C is asked
  again under the domain rules' definitions and conduct and the Outcome rule).
- **The tool descriptions:** a tool spec carries its advice apart from its description (`method`),
  sent with `'full'` and left out with `'off'`; the description keeps what the tool does and
  returns. On `dev-offcentre` eight pieces go: `find_gr_peaks`' "use first_shell.window_A for
  compute_bond_angles"; `list_sites`' "Use it to find the most disordered, anisotropic or anharmonic
  sites"; `compute_bond_angles`' "Use the whole first coordination shell (find_gr_peaks'
  first_shell.window_A): …" (its refusal stays stated on `window_confirmed`); `analyze_clusters`'
  "(take the window from find_gr_peaks: first_shell.window_A of that element pair)" and "Use it for
  cluster or polyhedron distortion … single-site displacement directions cannot show these";
  `analyze_modulations`' "and detect_symmetry" (cannot show a superstructure); `detect_symmetry`'
  "Compare tolerances to see which symmetry the local distortion breaks"; `measure_symmetry_lowering`'
  "which the displacement analyses cannot show — each site's cloud is centred on its own mean".
  With fit files, `compare_fit`'s "Use it for "where does the fit fail" …; run it once per dataset".
- **The checks:** a check whose text names an analysis, a tool or the window to use carries that
  advice apart (`methodCheck`), and arm C reads the fact and the verdict without it; a refusal still
  refuses. `compute_bond_angles`' cut-window refusal loses "Use the whole shell (find_gr_peaks'
  first_shell.window_A)." and keeps how the user confirms a window; its no-angles failure loses
  where the window comes from and keeps "Widen the window or check the elements"; the pooled
  directions note loses "Test single sites (site_ref), or analyze_clusters for a cluster's shape";
  `analyze_clusters`' window failures lose the `find_gr_peaks` window and its split-site note "Read
  the pattern …, the partial g(r) shell (find_gr_peaks) and the bond angles"; `find_gr_peaks`' two
  split-shell notes lose "Use first_shell.window_A, not the first peak alone" and "use it only if the
  user chooses it". Three reading hints that name no tool go too (decided by the manager on the
  user's delegation, 2026-10-19): `find_gr_peaks`' unseparated-shell warning loses "Ask whether the
  crystal has A–B distances there", `compute_bond_angles`' octahedron note "Read the shape from
  histogram_10deg_sin_corrected (1 = random directions), not from the raw modes", and
  `detect_symmetry`'s point-group note "Take inversion, chirality and polarity from this, not from
  memory". Two `analyze_modulations` notes say what a pattern is evidence of, an inference that arm
  C drops as it drops the domain rules' (decided by the manager on the user's delegation,
  2026-10-19; protocol 2026-10-19.19): the long-wave note loses "(common in RMC)" and the in-phase
  diffuse note "; as RMC's long waves and acoustic-like thermal motion give" (until protocol
  2026-10-24.2 "(as RMC's long waves …)"). Their facts and limits ("not a crystallographic
  superstructure"; the in-phase note's "(correlated motion), not order", until 2026-10-24.2 "not a
  superstructure") stay.
- **The evidence gate** runs and checks the same rules (decided by the manager on the user's
  delegation, 2026-10-19), but its send-back names no analysis and no tool: "… Your answer concludes
  "no order beyond thermal motion", but this conclusion needs its evidence from an analysis run in
  this question. Answer from the analysis's result (quote its numbers by reference). If none can
  run, say so in the answer instead of concluding." (arm B's names the analyses missing); where an
  analysis that would back it refused (protocol 2026-10-24.5), "… but an analysis that refused is no
  evidence: this conclusion needs its evidence from an analysis that answered in this question. …"
  (arm B's names the refusal and, for a window cut, both fixes). The record
  keeps the rule and the missing analyses for the reviewer (`gates`), and the end-of-turn finding
  names them too; the model reads neither (`harness/gate-send-back.arm-B.txt`, `.arm-C.txt`).

The other texts the model can receive name no analysis as method and are the same in both arms:
the Workbench's questions and their notes, the view rules (never offered in the live arms), the
loop's notes and reminders (`ask_user` is named there as the way to ask, not as a method),
`ask_user`'s own messages, and the "unknown tool" answer to a `read_skill` call in arm C (it lists the
tools offered). `modelVisible/harness.test.js` checks arm C's first
request for every skill name, the skills' step headings of more than one word, the prompt's method
sentences, `EVIDENCE_RULES`, the domain rules' inference parts and every piece of tool advice. So B against C measures what the
method adds: the skills together with the method in the prompt, the descriptions and the checks.

## Cases

Nine kinds, each a **dev** case (for improving tools and skills) and a held-out **test** case (for
reported numbers only): the same physics with other elements, directions, magnitudes or
structures, plus one dev-only case with fit data ([below](#the-fit-data-dev-case)) and the
judgment subset ([below](#the-judgment-subset-plan-item-16-part-ii)). Each case is
an `.rmc6f` supercell (reference sites and cell indices; coordinates wrapped into [0, 1) as
RMCProfile writes them) with its partial g(r). `run.js write` saves each as
`<out>/<case>/rmc_run/`, a run folder the Workbench opens, named as the harness names every run
(the Workbench names a run after its folder), with the truth beside it in `<out>/<case>/truth.json`.

The table gives each dev case. Its test case plants the same kind of order with other elements,
directions, magnitudes or structures; those values are sealed ([the sealed test
set](#the-sealed-test-set)) and appear nowhere else, this document included.

| Kind | Planted (dev) | The trap |
| --- | --- | --- |
| offcentre | Two-site disorder: Nb ±0.20 Å along [001] in KNbO₃ | A split site looks like large thermal motion in the ellipsoid alone |
| zone_boundary | Ordered antiphase Ti ±0.10 Å along c, q = (0,0,½); plus a stronger box-scale wave of the A cation | Identical single-site statistics to two-site disorder; the strongest mode is the decoy |
| breathing | Mo₃ triangles (O anion) whose size fluctuates collectively (sd 4 %), shape regular | Size is not shape |
| jahn_teller | Isolated CuF₆ (4 × 1.93 + 2 × 2.25 Å), long axis random over x, y, z | The first peak holds 4 neighbours; the average octahedron is regular |
| cell_wrap | Rigid V₂ dimers whose first site sits at fractional 0 (half its atoms written near 1) | The cell wrap |
| thermal_null | BaTiO₃ composition, thermal motion only (O flattened across its bond) | Literature off-centring; anisotropic but Gaussian O |
| symmetry_lowering | Every Ta shifted 0.07 Å along [001] (P4mm) | The 0.2 Å card still says Pm-3m |
| sites_only (achiral) | Zincblende CuBr, F-43m | No inversion does not make a structure chiral |
| sites_only (chiral) | B20 FeSi, P2₁3 | A genuinely chiral (Sohncke) group |

The paired dev cases: thermal rms 0.05–0.09 Å; supercells 6–8 cells a side, 1,024–2,560 atoms.

**The split is enforced in code** (`__tests__/split.test.js`, `__tests__/holdoutGate.test.js`):
- the test cases are sealed: only one loader may open them, on purpose
  ([the sealed test set](#the-sealed-test-set));
- the test set is pinned by a sha256 of its configurations, partials and truths;
- every paired dev kind has a test case of the same kind; a dev-only case (`devOnly`) outside the
  judgment subset has a kind no test case has;
- dev and test share no id, and no tool, skill, app source, doc or validation text names a test
  case's material, its planted values, its directions or its seed
  ([the holdout scan](#the-holdout-scan));
- `--split test` needs `--final` and a clean git tree, and every record names its commit;
- a supervisor finding on a test case may only report a bug, with the user's approval.

### The sealed test set

Since the holdout gate (2026-10-19) the test cases' definitions live in one sealed module,
`web_app/frontend/validation/holdout/testCases.js`, beside the test set's arm A / A+ / D records
(`holdout/arm-a-test.json`). Builders and agents never open, print, grep or import anything in
that folder. Before the gate the definitions sat in `cases.js` and this document stated their
values; git history keeps that text, and anyone who has read it is closed out of science,
threshold and scoring work on the validation.
- **The ids are public.** `cases.js` keeps `TEST_CASE_IDS` (ids only), so `splitOf`, the findings
  schema and the records can name a test case; nothing else of it is written outside the folder.
- **One loader.** `validation/holdoutLoader.js` `loadTestCases()` is the only module that imports
  the folder, by a dynamic import that runs only when asked on purpose: on `run.js`'s `--final`
  path (`arm-a`, `live` and `write` with `--split test`; the command line opens it with `--final`
  only, whatever the environment), or when the environment sets `RMC_VALIDATION_HOLDOUT=1`. Any
  other call throws `HoldoutSealedError` before importing. `casesFor(split)` is the async split
  helper; `cases.js` `casesOf('test')` throws.
- **CI opens it, locally it stays shut.** The frontend job's test step sets
  `RMC_VALIDATION_HOLDOUT=1`, so the tests that need the definitions run there: the test-set hash,
  the pairing with the dev set, the holdout scan's tokens, arm D on every case and the first
  requests of `caseBlind.test.js`. Without the flag they are skipped (`it.skipIf`), and the holdout
  scan looks for the public ids only.
- **The gate's own checks** (`holdoutGate.test.js`): no module of `src/**` or `validation/**` but
  the loader names the folder in code (a small lexer reads every string and template literal, so
  static and dynamic imports, `require`, globs and file paths are all caught, while comments may
  describe it); no test file and no model-visible snapshot writes a held-out id, and no test title
  is built from the held-out list; the loader and the command line refuse without `--final` or the
  flag.
- **Records written by the CLI.** A test-set file written inside the repository goes under the
  sealed folder (`run.js arm-a --split test --final --out validation/holdout/arm-a-test.json`) or
  the gitignored `validation/runs/` and `validation/cases-out/`; `run.js` refuses anywhere else.
- **Content-neutral.** The move changed no case: `buildCase` gives byte-identical output for every
  case, dev and test, before and after, and `TEST_SET_SHA256` holds. Nothing a run sends changes,
  so the protocol does not move.

### The fit-data dev case

**dev-fit-unconverged** (plan item 4.6, kind `unconverged_fit`, `devOnly` in `cases.js`) is the one
dev-only case outside the judgment subset. No other case has fit files, so without it `assess_convergence`,
`compare_fit` and the rmc-fit-assessment skill were never exercised. It has no test twin and
gives no reported number.
- **The configuration:** rock-salt MgO (a = 4.212 Å, 6 × 6 × 6 cells, 1,728 atoms), independent
  isotropic Gaussian displacements only (Mg 0.062, O 0.066 Å rms), seed 110. The structural
  answer is "no order", Fm-3m, centrosymmetric, achiral, non-polar.
- **The run folder:** besides `rmc_run.rmc6f` and `rmc_run_PDFpartials.csv`, RMCProfile's X-ray
  file family in the Demo's formats (`runFolder.js` `RUN_FILES`, `caseFiles`):
  `rmc_run_FQ1.csv` (Q, F(Q)_RMC, F(Q)_Expt), `rmc_run_FT_XFQ1.csv` (its Fourier transform: r(A),
  X_ray-calc, X_ray_exp_renorm) and `rmc_run-00.log` (400 rows, 10 minutes apart: clock time,
  moves accepted and generated, the χ² of F(Q)_1 and of X_ray_(R)1, Fortran E10.3). The
  `.rmc6f` header carries the run counters (216,000 moves generated, 16,257 accepted: 9.4
  accepted moves per atom). The run module detects the kinds (`xray_sq`, `xpdf`, `r_value`) as
  for a user's run. The "calculated" columns are the configuration's own X-ray G(r)
  (Faber–Ziman, f = Z) and its Lorch-windowed transform F(Q); the "experimental" ones carry the
  planted difference and noise (G(r) 0.004, F(Q) 0.002), from a second seeded stream (111).
- **The planted fit truth** (`truth.fit`, not scored):
  - *Not converged.* X_ray_(R)1's χ² falls fast, then slowly: 0.01 + 1.0·e^(−t/25) +
    0.02·e^(−t/600) per log row t. The plot's end looks flat, but ln χ² still drops 0.034 over
    the last tenth; the watchdog and `assess_convergence` say *improving*.
  - *One missed feature.* The measured G(r) has a Gaussian peak at 2.55 Å (height 0.15, sd
    0.08 Å, data above the model) in the empty gap between the Mg–O and the Mg–Mg/O–O shells,
    which the configuration does not reproduce. `compare_fit` puts 97 % of the G(r) misfit in
    r = 2.17–3.33 Å, its largest residual at 2.545 Å (model − data −0.15).
  - *One measurement in two spaces.* In F(Q) the same peak is a slow oscillation (period
    2π/2.55 ≈ 2.5 Å⁻¹) fading as 1/Q: 95 % of the F(Q) misfit lies in Q = 0.5–4.58 Å⁻¹. That is
    not a scale or box-size error, and the context marks the FT file as the same measurement.
- **How it is scored:** as every case, on its structure (no spurious claim and the four symmetry
  fields); `scoring.js` and the verdict schema are unchanged. The fit truth is exercised where
  the app checks it: the context (`convergence`, `datasets` with Rw, `partial_weights`,
  `configuration_optimization`), the two fit tools, and the answer check's convergence rule
  (an answer that says "converged" gets a finding and the Outcome "not verified", whatever its
  score). The record keeps every step and `truth.fit`, so a supervisor can read which fit tools
  ran and what the answer made of them. `__tests__/fitCase.test.js` pins all of this with a
  scripted model.

### The judgment subset (plan item 1.6, part ii)

On the paired cases the engines alone already score close to 1, so the case score cannot show a
model adding to them ([the bar after batch 1](#arm-a-results-no-model-2026-10-14-commit-b4fe09f)).
The judgment subset holds dev cases where A+'s fixed joining rules (`armA.js` `runArmSweep`) are
ambiguous or wrong, so a model that reads the tools with judgment has room to add.

- **The cases.** Each has `judgment: true` and a `situation` (`cases.js` `JUDGMENT_SITUATIONS`).
  They were `devOnly` until protocol 2026-10-19.11, when their held-out matches joined the test
  set (one per dev case, the same situation and kind, written blind from
  `validation/judgment-spec.md`; sealed since the holdout gate); `split.test.js` pairs them by
  situation and kind (in CI, where the loader opens the test set). The dev cases'
  parameters sit apart in `validation/judgmentCases.js`; `cases.js` appends them to `DEV_CASES` and
  builds them with its kinds. Every case before them builds byte for byte as before.

  | Situation | Case | Kind and option | Planted | The fixed rule at issue |
  | --- | --- | --- | --- | --- |
  | Two features at once | dev-jt-off-centre | `jahn_teller`, `offCentre` | Jahn–Teller octahedra whose centre is also off-centred: both on the centre | A split centre explains its own split shell |
  | A signal at 2.5–4σ | dev-weak-order-a, dev-weak-order-b | `zone_boundary`, `weak: true`, `decoy: null` | A weak cell-doubling order of B, below half its thermal rms; the truth names q only. One seed on each side of the bar | The fixed 3σ bar |
  | Two tools read differently | dev-eight-site | `offcentre`, `wells: 8` | Eight-site <111> off-centring of B | A split only from the ellipsoid (a resolved PC1, `isSplitSite`); the 3 + 3 shell it leaves unexplained |
  | A choice only the user can make | dev-octahedral-clusters | `elongated_clusters` (new) | M6 octahedra, each elongated along a random axis: a cluster distortion | None: the control (below) |

  `elongated_clusters` is the subset's one new kind: no other kind plants a cluster distortion. The
  other cases are options of existing kinds; without the option a kind draws exactly as before.
- **The control.** The octahedral clusters' M–M shell is split across a trough. The Workbench joins
  it (an octahedron) and asks when an analysis' window cuts it. The choice of window is the user's,
  and the case's truth follows the default, the whole shell. Since protocol 2026-10-19.8, A+ takes
  its cluster windows from the same shell reading (`armA.js`), so its fixed rule and the default
  coincide: on the score this situation is a control, and the judgment (saying which window was
  used and what the other gives) shows in the answer only. The case was designed against the
  earlier rule, a cut at the first trough. It still pins the question and its default
  (`judgment.test.js`) and would show a reader that leaves the one shell reading again.
- **Seeds.** The first three cases take the first seed, upward from their base, on which A+ errs
  only where the situation says, so each measures its own situation. The situation's own failure
  held on every seed tried (the octahedral case's under the cluster window before protocol
  2026-10-19.8). The side failures seen are reported with the PR as findings for the engines. The
  weak-order pair shares one amplitude, set so the expected significance is 3σ: since protocol
  2026-10-19.19, the median realised σ against the modulation engine's independent-cell null, by the
  spec's formula and checked with the engine on fresh configurations. Its seeds come from an
  outcome-neutral rule: one realised on each side of the bar (`judgmentCases.js`).
- **Scoring is unchanged.** Three things to read with it:
  - the verdict has no confidence field, so an honest "marginal" shows in the prose, the answer
    checks and the computed verdict, not in the score;
  - stating the default of a window the user skipped is not scored; the planted truth is the one
    the Workbench's default (the whole first shell) leads to;
  - each case keeps its kind's tolerated partial descriptions, and the new kind tolerates its own.
- **Dev only, until the held-out version.** A fresh author writes the held-out matches from
  [`validation/judgment-spec.md`](../../web_app/frontend/validation/judgment-spec.md) alone: the
  situations, the kinds and options, the knobs, and no value. They do not read `judgmentCases.js`,
  the dev records or the judgment pins. The lead checks their cases against the dev ones without
  showing them. The matches join the test set as a new pinned version on the user's decision
  (done at protocol 2026-10-19.11; their values are sealed with the other test cases). Until
  then:
  - no test case is a judgment case, and `TEST_SET_SHA256` is unchanged (`split.test.js`);
  - a judgment case never removes a test case's token from the holdout scan (`holdoutScan.js`);
  - `__tests__/judgment.test.js` scans the spec, this section and the change note for the dev
    cases' ids and values.
- **How the report reads it.** Each split with judgment cases gets a table: each column's mean
  score per situation, over the subset and over the rest, every run counted (`report.js`
  `renderJudgment`; a record is a judgment case's when its truth says so, `truth.judgment`). The
  comparison beside the case score (part i) lists their components like any case's, and a
  transcript's planted truth names the situation.
- **Two situations through the harness** (`judgment.test.js`, a scripted model): a window that
  cuts the octahedral clusters' split shell is asked about, and the simulated user's skip keeps
  the whole shell; a superstructure claimed below the bar scores as a claim while the answer
  check flags it and the computed verdict is "not verified".
- **Where the results are.** In the PR, `validation/baseline/arm-a-dev.json` and the `arms.test.js`
  pins, not here: the held-out author may read this document.

## The holdout scan

`__tests__/split.test.js` (with `__tests__/holdoutScan.js` and its lists, `__tests__/holdoutLists.js`)
keeps the test cases out of what the tools, the skills, the docs and the validation's own text
say. It derives each test case's tokens at run time from its parameters, read through the loader
([the sealed test set](#the-sealed-test-set)), so no test-case value is written in the test; a
failure prints the line with what matched masked (■). Every kind of token needs the parameters,
so the full scan runs where the loader opens them (CI's `RMC_VALIDATION_HOLDOUT=1`); without the
flag it looks for the public ids only.

- **Tokens, from the parameters:**
  - the id, anywhere;
  - the composition (two or more elements): a formula with exactly that element set, with or
    without subscripts, or, for three or more elements, an element-dash chain naming each once;
  - the cluster written as a formula of its own (a tetrahedron's M₄, a dimer's M₂);
  - a structure-type word that only a test case uses;
  - each planted value within six words of its element: every length in Å (displacement, bond
    lengths with their split and mean, cluster edge and elongation, the dimers' shared and own
    shifts, the lattice parameter), a thermal width the case sets itself, and a breathing sd (as a
    fraction and in per cent); a lattice parameter and a thermal width count next to a cation of
    the case, not next to its anion (the perovskite cases share O);
  - the supercell (N³, N × N × N, [N, N, N]) within six words of its element;
  - each planted direction family ([uvw] or <uvw>, any order and sign) and wavevector
    (q ≡ −q ≡ q + G, any order) within six words of its element;
  - the seed, a whole number within six words of the word "seed".

  Values, supercells, directions and wavevectors count only where they differ from the dev cases
  of the same kind: those are the kind's physics, not the test case. A seed always counts. A
  thermal width the kind fixes for every case (the perovskites') is no case's parameter and never a
  token. A judgment dev case ([above](#the-judgment-subset-plan-item-16-part-ii)) is never that
  comparison: it cannot remove a token. `split.test.js` checks that every test case yields its seed
  and a token for each value and supercell its dev twins do not share.
- **Not hits:** a bare element symbol (they are everywhere); a formula with more or fewer
  elements; a chain that repeats an element (A–X–A is a bond angle, which every study of that pair
  writes); a number with more digits (0.371 is not 0.37) or more than six words from its element;
  a supercell next to the word "supercell" alone (a real run's, the Demo's among them, is not a
  case).
- **Where:** `web_app/frontend/src` (js, jsx, md, and the snapshots of what the model sees), every
  file in `.claude/skills`, `rmc_toolkits` (py), every md under `docs/` (this document, the change
  notes in `docs/changes/unreleased/`, the CHANGELOG and [AGENT_PLAN.md](../AGENT_PLAN.md) among
  them), the top-level md files, and the validation itself outside the sealed folder: its js and md
  files (`validation/**`, tests and `judgment-spec.md` included) and its findings (the decision
  records and the batches). The sealed folder is never read. A test asserts this coverage.
- **Checks on the scanner:** every test case yields more than its id, and the scan finds every
  test case in its own truth (summary and traps, prose written from the same parameters). Its own
  rules are pinned on made-up cases (`holdoutScan.test.js`).

Every hit must sit in one of two lists in `holdoutLists.js`, and each entry pins its exact number
of hits, so a new mention, or a removed one, fails until the entry is reviewed. Any hit in neither
list fails the test. An id-only entry (`case: '*'`) covers every held-out id in its file together;
it allows ids and nothing else.

**Open leaks, pending the user** (`HOLDOUT_LEAKS_PENDING_USER`). Text the model can read in arm B
(the system prompt, the tool descriptions and results, `read_skill`'s output, every skill file)
that names a held-out test case. Such a hit is a leak, not a known risk: the test prints it on
every run and pins it line by line, and CI stays green only so that the decision stays with the
user, whose text the skills are.

**The list is empty since protocol 2026-10-19.2.** It held three lines of
`.claude/skills/coordination-shells/lacunar-spinel-shells.md` (`:9`, `:31`, `:32`): the test
breathing case's cluster was the cluster of the skills' reference compound, whose published
distortion that reference file describes. The test case was recomposed instead of the skill
(Protocol 2026-10-19.2, below), and nothing names its new composition: not the skills on main,
not PR #71's drafts (read-only scan, 2026-10-19), not the snapshots. The list and its checks stay
for any future skill-text hit.

**Known risks: the allow-list** (`HOLDOUT_ALLOWED`). Only text no model can read, in any arm: the
test checks that no allowed file is under `.claude/skills`, among the files `copilotSkills.js`
bundles for the prompt and `read_skill`, or in the snapshots of what the model sees; that no app
module (`src` outside `__tests__`) names its path, so none can import it; and that no snapshot
contains its matched text.

| File | Case · token | Hits | Why it is not model-visible |
| --- | --- | --- | --- |
| `web_app/frontend/src/__tests__/fixtures/symmetryStructures.js` | test-symmetry-lowering · composition | 1 | Test only: a symmetry-finder fixture of a real compound, imported only by the symmetry tests; unrelated to the planted physics |
| `docs/algorithms/copilot-validation.md` | every test case · id | 10 | This document names the cases by id; ids only |
| `web_app/frontend/validation/cases.js` | every test case · id | 15 | `TEST_CASE_IDS`, the public manifest of the held-out ids (one id also sits inside another) |
| `web_app/frontend/validation/findings/decisions/2026-10-19-test-breathing-recomposed.json` | every test case · id | 2 | A decision record on a test case names its case |
| `web_app/frontend/validation/__tests__/holdoutLists.js` | every test case · id | 2 | The lists themselves: each entry names the case and the file it covers |

The first entry was found on 2026-10-19 and reviewed against that rule; it points at a real
compound, so it remains a known residual of the seal. Four entries found that day for the test
breathing case's old cluster (the llm README, a triplets engine test, ai-assistant.md and
bond-geometry.md, all doc or test only) left the list with protocol 2026-10-19.2: they name the
skills' reference compound, which is no longer a test case. The id-only entries came with the
holdout gate, when the scan first read this document and the validation.

## The question and the verdict

Every case, arm and model gets the same open question (`CASE_QUESTION`): is there order beyond
the average structure, and what is the space group (chiral, centrosymmetric, polar)? The live arms
end their answer with one `VERDICT {…}` line in a fixed schema (`verdict.js`): features of seven
kinds (`off_centre`, `modulation`, `cluster_distortion`, `cluster_breathing`, `rigid_cluster`,
`jahn_teller`, `symmetry_lowering`), each with element, direction [uvw], magnitude and q, plus the
four symmetry fields. Arm A writes the same object from its readings. Nothing is parsed from
prose: an answer with no parseable verdict scores 0. Since protocol 2026-10-17.2 a number in the
verdict may be a reference to a tool result (`"magnitude_A": {{r3.sites.0.shift_A}}`): the verdict
is parsed from the answer with its references filled, and one that names nothing leaves the JSON
unparseable (score 0), as an invented number would be wrong.

## Scoring (`scoring.js`), fixed before any run

The case score is the mean of these components, each 0 to 1:
- each planted feature: 1 if found (kind and element) and accurate, 0.5 if found with a quantity
  off, 0 if missed;
- "no spurious claim": every claim is planted or tolerated (a true partial description, e.g.
  "F off-centre" in the JT case); in a no-order case this is the **over-claim** test;
- each symmetry field the truth defines: wrong or unanswered scores 0.

Tolerances, from the physics of the cases:
- **Direction 15°:** half the smallest angle between distinct low-index families (<110>–<111> is
  35.3°). ±u count as one for a split site; any cubic image counts for a JT axis.
- **Magnitude: 25 % or 0.03 Å, whichever is larger.** 0.03 Å is half the smallest thermal rms of
  the cases, below which a static shift cannot be told from the thermal width. 25 % accepts
  √(rms1² − rms⊥²) and fails quoting the thermal-broadened rms.
- **Wavevector 0.02 r.l.u.:** with 8 cells a side the allowed q are multiples of 0.125; q ≡ −q ≡ q + G.

Reported beside the score, per run:
- the answer checks against the contradictions the score found (caught, missed, not covered);
- unsafe actions (a window marked as the user's choice when no question reached the user; must
  be 0): a `compute_bond_angles` step with `window_confirmed: true` and no earlier `ask_user`
  step, the model's or the Workbench's, whose question reached the user (the step's `choice` is
  set: an option, a typed answer or a skip). A malformed or refused question (`result_status`
  `error` or `skipped`) does not count, nor a window the Workbench settled with no card, which
  records no `ask_user` step (since protocol 2026-10-19.13; before it any earlier `ask_user` step
  counted, whatever its result);
- unknown-tool, failed and refused calls, read from each step's `result_status`: unknown-tool
  calls end `unknown_tool`; failed calls are every error the model read instead of a result
  (`error`, `unknown_tool` and `skipped`); refused calls (`refused_calls`, added with the tool
  transcript) are the tools' own refusals (a `fail` check, no error), which count as analyses that
  ran;
- **number provenance:** each number with a decimal point or a unit, matched to a value the
  model was told, or "unsourced": per step, the tool message it read (the data as fitted, the
  checks, the error, then the note; the first that holds the number), the run context as it read
  it and the question. A Workbench question sends no message: its note, the text the model read in
  its place, counts, never its card's data or checks (since protocol 2026-10-19.13; before it every
  step's data and checks counted, the card's included, and no note or error);
- **reproducibility over 3 runs:** an identical verdict signature (kinds, elements, q, group,
  symmetry facts), pair agreement, and the score's sd;
- **the evidence gate** (THY-18): how often an answer was sent back because a conclusion's
  analysis had not run, and which runs still lacked evidence at the end;
- **numbers by reference:** references written, filled and unknown (provenance traces the filled
  answer, so a number written by reference always has its source);
- **the computed verdict against the score:** *earned* (achieved, score 1), *over* (achieved,
  score below 1: what the checks cannot see), *under* (not achieved, score 1) and *held* (not
  achieved, score below 1: the verdict warned). The Outcome is the app's, computed from its checks;
  the model only explains it.

Per split, the report also takes the case score apart: its parts (the feature-only score among
them), every component against A+, and every disagreement with the planted truth
([below](#the-report-beside-the-case-score-plan-item-16-part-i)).

**Protocol 2026-10-17.2 (THY-18).** The live arms run with the evidence gate, numbers by reference
and the computed verdict, exactly as the app does. The case score is unchanged; the three are
reported beside it (`report.js`: one more table per split) and in each transcript (the answers the
gate sent back, the answer as written and as filled, the computed outcome with its reasons). Since
architecture candidate #4 (2026-10-19) the three are one module, the verification; the harness
adds the verdict schema's words to the claims its gate reads (`VERDICT_PHRASINGS`), which the
app's catalogue read itself before. What a run sends and records is unchanged (no protocol bump).

**Protocol 2026-10-19.1 (the case-id leak).** Until this version the harness named each run, and
its `.rmc6f` and partials files, after the case id. The run name reaches the model as the
context's `run`, in the first message, so every arm B and C run was told its case kind
(`"run": "dev-zone-boundary"`). Runs are now named `rmc_run`, the same for every case. Live runs
under an earlier protocol, if any, are not comparable with later ones and must be run again. Arms
A, A+ and D have no model and are unaffected; their records keep protocol 2026-10-17.1.

**Protocol 2026-10-19.2 (test-breathing recomposed).** The held-out breathing case planted the
cluster and anion of the skills' reference compound, whose published shells
`.claude/skills/coordination-shells/lacunar-spinel-shells.md` gives and `read_skill` serves to arm
B. That leaked the case's kind to arm B only. The fix goes in the test case; the skill is the
user's text and is unchanged. Since the holdout gate the new case's values are sealed with the
other test cases ([the sealed test set](#the-sealed-test-set)); see `TEST_SET_SHA256`.
- **The new composition:** a cluster element and an anion a fresh author chose without seeing dev
  results, dev parameters or arm records, with a lattice parameter and an edge that keep the old
  shell separation; shape, breathing sd, thermal sd and supercell unchanged. Neither element, the
  cluster's formula, the pair's chain nor the elements' names appears in any skill (main or PR
  #71's drafts) or in the model-visible snapshots.
- **The seed, by a rule fixed before looking:** the first seed of a fixed sequence whose drawn
  breathing sample sd over the case's clusters (the generator's own `cellState` draws) lies in a
  fixed window about the planted sd, at most 20 tries; the second seed tried was taken, and no
  other case uses it. The rule depends only on the case's own draws, not on any model or arm
  result.
- **Difficulty, honestly:** the planted and the realised radius spreads stay within a few per cent
  of the old case's, and A+'s cluster engine reads the breathing more strongly than before (the
  numbers are sealed with the case).
- **Approved by the manager on the user's delegation, 2026-10-19**
  (`validation/findings/decisions/2026-10-19-test-breathing-recomposed.json`). The user may reverse
  it: keep their compound in the test case and disclose it as a known arm-B advantage.
- **Nothing is lost:** no arm B or C run has touched the test set. `TEST_SET_SHA256` changes, so
  any result on the old test set would not be comparable.
- **Arms A and A+ on the test set**, re-recorded with `--final` at 4ab9eab: only the three
  test-breathing records change (readings, the verdict's element, the truth), and no score moves.
  A scores 0.50 (no page sees cluster breathing). A+ scores 1: the breathing found, no shape
  distortion. D scores 1. The means stay A 0.81 and A+ 0.99. The records are in
  `validation/holdout/arm-a-test.json` (sealed). The dev records are unchanged: a re-run equals
  `arm-a-dev.json` record for record, apart from protocol and commit.
- **The dev set overlaps the skills, and may.** dev-breathing (Mo₃ triangles, O anion) was checked
  against the same rule. No Mo₃, Mo–O or Mo appears in the skills on main or in the snapshots.
  PR #71's draft `local-order/order-signatures.md:25` reads "Breathing Mo₃ triangles (size sd
  4 %)": dev-breathing's cluster and its planted breathing. The dev set is for development, and
  the skills learn from it. Once that draft merges, read an arm B result on dev-breathing with this
  in mind; the reported numbers come from the test set.
- **The holdout scan's open-leak list is now empty**, and the four allow-list entries for the old
  cluster have left it ([the holdout scan](#the-holdout-scan)).

**Protocol 2026-10-19.3 (`list_sites` shows the split, plan 4.1).** Each site row gains its
per-axis excess kurtosis on resolved axes (`kurtosis_axes`), PC1 as a direct-lattice [u v w]
(`pc1_uvw`), and, for a site the tools call split (one rule, `isSplitSite`), its static
displacement `delta_A` with its sampling error `delta_se_A`; the tool's description says so. What
the model sees changes (the `list_sites` result and its description), so live runs under earlier
protocols are not comparable. Arms A and A+ do not call `list_sites` and are unchanged. Approved by
the manager on the user's delegation, 2026-10-19 (`validation/findings/decisions/`).

**Protocol 2026-10-19.4 (the turn module).** The harness had assembled each question itself,
beside the dock's own assembly, and the two had drifted apart. Both now run the turn module
(above), and four things change for the harness; each is pinned in `harness.test.js`:
1. **The warmer retry.** After a reply that kept repeating itself, the loop asks again warmer: the
   dock sent that request at a temperature of 0.6 or more, the harness at the model's own (0.2), as
   `streamCopilotChat` ignored the loop's flag. It now runs at 0.6 or more in both
   (`providerModel`, `ANTI_REPEAT_TEMPERATURE`).
2. **The no-tools prompt.** A model that refuses tools gets the tool-free system prompt, as in the
   dock; the harness had kept the tool rules.
3. **Provenance against what the model read.** Each answer number is matched against the run
   context as the model was shown it (trimmed to its budget, the turn's `start.shown`), no longer
   the whole context: a value the trim removed was never shown. No case's context is over the
   budget at the default view (the largest is ~4,500 of 6,000 characters, ~5,600 for the Demo with
   the PCA table), so no default record moves. (Since protocol 2026-10-24.1 the budget is 6,500
   characters and the convergence block carries `last_10pct`: the fit case ~4,200, the Demo ~6,100.)
4. **Every request as sent.** The record keeps `requests` (below).

Changes 1 and 2 change what a run sends, so the protocol moves. Live runs under an earlier
protocol that looped or met a model without tools are not comparable; arms A, A+ and D have no
model and are unaffected.

**Protocol 2026-10-19.5 (the run policy, plan item 1.2).** What the model is sent is unchanged:
the model-visible snapshots did not move. Two things change, both [below](#the-run-policy-plan-item-12):
1. **An LM Studio run starts only at 48k or more.** The pre-run check reads the loaded window
   and refuses a run below 49,152 tokens. A turn that a 32k window would refuse, a score of 0 for
   the harness's window rather than the model, can no longer happen in a run.
2. **Each record states its run policy** and the window it ran under (`record.policy`, and each
   request's `num_ctx` on Ollama).

The protocol moves because the run's conditions are now part of it: a record under an earlier
protocol does not say which window it ran under, so it cannot be set beside one under this. No arm
B or C run exists under an earlier protocol; arms A, A+ and D have no model and are unaffected.

**Protocol 2026-10-19.6 (plan item 4.6: the fit-data dev case).** dev-fit-unconverged joins the dev
set ([Cases](#cases)). The question, the scoring, the verdict schema, the tools and the skills are
unchanged; the run module reads the new files exactly as it reads a user's run.
- **Arms A and A+ on it**, appended to `arm-a-dev.json` (`arm-a --split dev` at 2383578): both
  score 1. Neither claims anything, and both read Fm-3m with its facts, A from the 0.2 Å card and
  A+ at the noise floor (0.034 Å). A+'s cluster engine finds extended Mg and O networks, not
  finite clusters, and its modulation test no zone-boundary order. D scores 1. The re-run's other
  27 dev records equal the committed ones apart from protocol and commit, and stay as recorded.
- **The means:** over the nine paired dev kinds they stay A 0.81 and A+ 0.99; over all ten dev
  cases they are A 0.825 and A+ 0.99. The test set and `TEST_SET_SHA256` are unchanged.
- **Arm B on dev:** the new case's first request is the dev set's largest (11,178 tokens: the two
  fit tools and four more context sections). `contextBudget.test.js` now pins its worst turn:
  35,906 tokens needed, 3,138 over 32,768 (dev-symmetry-lowering's was 35,119, 2,351 over). The
  snapshots of what the model sees (the Demo, dev-offcentre) are unchanged.
- **What the model is told:** nothing about the case. Its id, its kind and the id's stem reach no
  request and no tool result (`caseBlind.test.js`, which now also checks that the dev cases
  between them offer every tool its sweep calls: the fit tools only here).

**Protocol 2026-10-19.7 (the model port, architecture candidate 3).** The provider quirks moved
behind the model port (`providerModel`, [the llm README](../../web_app/frontend/src/llm/README.md),
*The model port*). Three things change for a run, each pinned (`runPolicy.test.js`,
`modelPort.test.js`, `turn.test.js`):
1. **An Ollama run's final request carries no tools.** Ollama's native `/api/chat` has no
   `tool_choice`, so `'none'` was not sent and the tools were still offered. A final reply of calls
   only then ended on the loop's "budget used" answer, scored as no verdict. The port now enforces
   `'none'` by sending that request no tools; the earlier calls and results still travel as
   messages. **This is what an Ollama run sends in its last round, so Ollama records under an
   earlier protocol are not comparable on their final round.** The OpenAI-compatible servers
   (LM Studio, the cloud) get the same body as before.
2. **Each request as sent keeps the context guard's plan** (`sent.plan`, below): the server kind,
   Ollama's `numCtx` and `numPredict`, LM Studio's loaded context, the estimate, and the tool
   choice's form. Each request's `num_ctx` is read from it.
3. **A request the guard refuses is not recorded as sent.** The guard decides in the port's
   `prepare`, before anything is sent. A refused request has no entry in `requests`, as one the
   spend cap refuses, and the run's `error` is the `ContextLimitError` with its sizes. It had been
   recorded as sent, followed by the error. The choice keeps `requests` meaning "what reached the
   model": the refusal is the run's outcome, not a request, and its numbers are in its message.

What the model sees in every first request is unchanged (the model-visible snapshots pass as they
were). No arm B or C run exists under an earlier protocol; arms A, A+ and D have no model and are
unaffected.

**Engine verdicts (2026-10-19, architecture candidate 5; not a protocol change).** The cluster and
modulation engines decide what their statistics show and return it (`verdict`; ai-assistant.md
Steps 18–19); `analyze_clusters` and `analyze_modulations` only render it, and what the model sees
is byte-identical (the model-visible snapshots), so live runs under protocol 2026-10-19.7 stay
comparable. A+ now calls the engines (`run.engine('clusters' | 'modulations')`) and maps their
answers: distorted → `cluster_distortion`, breathing → `cluster_breathing`, rigid →
`rigid_cluster`, a cell doubling or order at a general q → `modulation` at the q that decided it.
A+ now reads clusters by the tool's rules, where it had its own: trigonality counts as a
distortion of tetrahedra; rigid is blocked by the engine's split member sites (at least half of
the cluster's sites split by `isSplitSite`) instead of arm A's split site (a prolate PC1 with
κ₁ ≤ −0.5 on any site of the element), and by a pattern (an axis family, a trigonal lean) beside a
narrower edge spread; collective breathing no longer blocks it (the edge spread is read against the
shape null, which keeps each cluster's breathing). On the dev set every score and every claim is
unchanged (the `arms.test.js` pins hold); only the ten A+ records' cluster and modulation readings
change, from the tool's sentences to the engines' answers; they are re-recorded in
`arm-a-dev.json` (at 3f90a80; the A and D records stay as recorded). Only a held-out case would
exercise the rule changes; the test set was not run, and its A+ records are re-recorded at the
freeze. Approved
by the manager on the user's delegation, 2026-10-19 (`validation/findings/decisions/2026-10-19-engine-verdicts.json`).

**Protocol 2026-10-19.8 (one shell reading per pair, per run, architecture candidate 6).** The first
coordination shell of every pair is read once (`src/llm/context/shells.js`, ai-assistant.md Step 9),
and the context, `find_gr_peaks`, the Workbench's questions, the evidence gate, the answer checks,
`compute_bond_angles`' refusal and A+ all read it; before, five readers searched their own ranges and
integrated their own way. What the model sees changes:
1. **The first peak alone ends at the minimum between the first two resolved peaks** (approved by
   the manager on the user's delegation, 2026-10-19), no longer at r ± FWHM, whose half-height walk
   reached into the second peak: the Demo's Ta–Se first peak is 2.21–2.65 Å holding 3.25 Se per Ta
   (2.33–2.71 Å and 3.86–3.9 before; the known answer, a 3 + 3 split with overlapping widths, gives
   3 ± 0.1 where r ± FWHM gives 3.4–4.0). A shell of one peak has none.
2. **`find_gr_peaks`:** `first_peak_window_A` leaves; the first peak is `first_shell.first_peak_only`
   (window, CN, trough) for every shell of two or more resolved peaks, and the split-shell note
   gives its CN; the shell no longer depends on `r_max` (a pair whose first peak lies above it gets
   a note and the shell, where it got "no peak"), and the tool's description says so; the
   constraint warning loses its FWHM clause (nothing moved on the Demo or the dev cases).
3. **The context's `pair_correlations`:** peaks below 8 Å (the Demo's Ga–Ga 7.33 Å, Ga–Ta's second
   peak), `first_peak_to_A` / `first_peak_cn` for every shell of two or more resolved peaks (the
   Demo's Ta–Se and Se–Se, dev-offcentre's K–Nb and Nb–O), and each CN integrated over exactly the
   window it quotes (third-digit moves: the Demo's Ta–Ta 2.76 → 2.75, Ga–Ta 12 → 11.9,
   dev-offcentre's O–O 8 → 7.99).
4. **The Workbench's "First peak only"** on the Demo: 2.21–2.65 Å, CN 3.2 (was 2.33–2.71 Å, 3.9).

The snapshots moved where expected: `find_gr_peaks.json`, `dock/first-request.txt` and
`harness/first-request.arm-B.txt` (the context), and the two `tools` files (the description); the
`compute_bond_angles` snapshots did not (the cut rule is unchanged; only its no-angles advice now
points at `first_shell.window_A`). The dock's worst turn grows by 90 tokens (the context +73, the
description +16: 38,652 needed, 10,500 to spare at 48k; `contextBudget.test.js`), arm B's worst
(dev-fit-unconverged) by 16 (35,922 needed, 13,230 to spare). Live runs under earlier protocols are not
comparable on shell, coordination or bond-angle questions. **Arms A and A+ on dev:** every score
and claim unchanged; three A+ records move, only in the window their cluster reading names, whose
low end is now the catalogue's (onset − one grid step, 0.01 Å) — dev-breathing Mo 2.40 → 2.39,
dev-jahn-teller F 2.50 → 2.49, dev-symmetry-lowering O 2.52 → 2.51 Å — with the engines' verdicts
unchanged; re-recorded in `arm-a-dev.json` at f763969, the other 27 records as recorded. The test
set was not run.

**The skill review is a precondition for the freeze.** The skill text is the user's and is not
edited: until the user updates the coordination-shells skill, arm B reads (`read_skill`) "the
first-peak window 2.33–2.71 Å holds 3.9 Se per Ta" (SKILL.md step 2, and the Demo table in
`lacunar-spinel-shells.md`) while `find_gr_peaks` says 2.21–2.65 Å and 3.25. The skill lint's new
value check finds the four places and holds them on its allow-list with both values.

**The cut rule did not change; a recommendation waits for the user.** Two readings of a window that
cuts a shell live in the shell module: the Workbench asks on the CN margin (the window holds fewer
atoms than the shell by ≥ max(0.5, 15 %)), `compute_bond_angles` refuses on the counts (more than half
the central atoms unbonded, or fewer than 60 % with the most common count). `run.js shell-cuts`
tabulates both on every natural window of the dev cases and the Demo (the whole shell, its first
peak alone, its top 0.2 Å earlier and later, and the Demo's historical Ta–Se windows). Of 149
windows they agree on 134 (6 cut by both, 128 by neither):
- **The counts alone refuse 11:** dev-offcentre's whole Nb–Nb shell (3.40–4.24 Å, CN 5.46; 52 %
  share the modal 5 — its counts spread because the shell runs into the next, 4.2 % at the
  boundary), and 10 windows with the top 0.2 Å early whose CN stays within the margin (e.g. the
  Demo's Ga–Ga 6.88–7.53 Å, 11.2 of 11.9; dev-jahn-teller's Cu–K 6.87 of 8 at 38 %).
- **The CN margin alone asks on 4:** first peaks whose counts form a sub-shell — the Demo's Ta–Se
  (3.25 of 6; 72 % have exactly 3), dev-jahn-teller's four short Cu–F (3.98 of 6; 92 % have 4), the
  Demo's Se–Se (2.78 of 12; 76 %) and dev-offcentre's one short Nb–O (0.89 of 6; 72 %).
- **Both:** the spot check's 2.3–2.7 Å and the old first-peak window 2.33–2.71 Å on the Demo, two
  first peaks and two early cuts on the dev cases.

The recommendation (decision record `validation/findings/decisions/2026-10-19-shell-catalogue.json`,
not approved): `compute_bond_angles` refuses on the CN-margin reading, the Workbench's own, so a
window the Workbench did not question is not refused and a whole shell never is; the 10 slightly
cut windows would then run (their CN within 15 % of the shell's), and the 4 sub-shell first peaks be
refused unless the user chose them — which the Workbench's question already asks. The user decides.
*(Approved by the user and implemented on 2026-10-19, with exactly those 15 changes: protocol
2026-10-19.18 below, `validation/findings/decisions/2026-10-19-shell-cut-rule.json`.)*

**Protocol 2026-10-19.9 (plan item 1.6, part ii: the judgment subset).** Five dev cases join the dev
set ([the judgment subset](#the-judgment-subset-plan-item-16-part-ii)). The question, the scoring,
the verdict schema, the tools and the skills are unchanged, and the run module reads the new cases
as any case. What a dev batch holds changes, so the protocol moves: a dev batch under an earlier
protocol has no judgment case, and its means are not comparable with one under this.
- **Arms A, A+ and D on the new cases** are appended to `arm-a-dev.json`, under this protocol and
  on top of protocol 2026-10-19.8's re-recorded windows. The earlier dev records stay as recorded:
  a re-run equals each of them apart from protocol and commit. D scores 1 on each new case.
- **The test set is unchanged:** every test case builds byte for byte as before, and
  `TEST_SET_SHA256` holds.
- **What the model is told:** nothing about the case. No new id, kind or id stem reaches a request
  or a tool result (`caseBlind.test.js` covers every case). The model-visible snapshots are
  unchanged. Arm B's largest first request on dev is still dev-fit-unconverged's
  (`contextBudget.test.js`).
- Approved by the manager on the user's delegation, 2026-10-19
  (`validation/findings/decisions/2026-10-19-judgment-dev-cases.json`).

**A+ reads the app's rules (2026-10-19, architecture candidate 8; not a protocol change).** A+ had
three rules of its own where the app has one. It now calls the app's functions, which
`SWEEP_RULES` holds, and a test checks they are the app's:

| Rule | A+ before | The app's, which A+ now applies | Evidence |
| --- | --- | --- | --- |
| Split site | arm A's: PC1/PC2 ≥ 1.3, PC1 resolved, κ₁ ≤ −0.5; δ = √(rms1² − rms2²) | `isSplitSite`: PC1 resolved, κ₁ ≤ −0.5 and at least 3 sampling errors (√(24/n)) below zero; δ ± its error from `list_sites` (`splitDisplacement`, √(rms1² − rms⊥²)) | On synthetic prolate Gaussian sites arm A's rule called 3.6 % (216 copies) and 14 % (100 copies) split, `isSplitSite` none (2,000 seeds each). A planted ±0.15 Å split over a 0.05 Å width at 216 copies is found in all 300 seeds, its δ within 1.96 errors in 94 % of them. rms2 alone biases δ low (ai-assistant.md Step 20). |
| Cluster elements | every same-element first shell peaking below 3.6 Å | the evidence gate's: the element's own first shell holds 0.5–4.5 neighbours (`clusterElementsOf`, `CLUSTER_CN_MAX`, now in `shells.js`; `evidence.js` imports it) | Dimers 1, trimers 2, tetrahedra 3, octahedra 4 are finite clusters; CN 6–12 is a network. On the ten paired dev cases the CN rule tries the planted V₂, Mo₃ and the F around Cu. The eight elements the old rule also tried there (the O of four perovskites, FeSi's Fe and Si, rock-salt Mg and O) are extended networks that the engine itself refuses. |
| Lower group's tolerance | the noise floor | `lowerGroupTolerance`: `NOISE_FLOOR_MARGIN` (1.25) × the floor, `measure_symmetry_lowering`'s default; the Workbench's tolerance question counts a group only where it holds to the same margin (`clarify.js`, until now its own `NOISE_MARGIN`) | No dev group moves: every lower rung ends below 0.7 × the floor and every true group holds past 1.25 × it. The cost: a lowering whose group holds only to 1–1.25 × the floor is not read. A synthetic shift of 0.0147 Å, whose group holds to 1.18 × the floor, now reads Pm-3m, while one of 0.027 Å (2.2 ×) reads P4mm. The app chose the margin, and no dev case shows it wrong. |

Arm A's rules are the frozen protocol's and are unchanged. What A+ keeps of its own is only how it
joins the engines:
- one claim per element;
- an ordered mode explains its element's split sites;
- an off-centred centre explains its own split shell;
- a mode's polarisation was the PC1 of the element's most prolate site (PC1/PC2 ≥ 1.3), when the
  app had no polarisation rule (plan section D); since protocol 2026-10-19.12 A+ reads the
  modulation engine's polarisation (below).

**Arms A and A+ on dev** (`arm-a --split dev`: before at a2c5dad and after at 63569f3 on the ten
paired cases; after again at 3de8fc5 with the judgment subset):
- Every arm A and D record is identical, and so is every score and claim.
- Every A+ record moves in its readings. On the ten paired cases:
  - The group is read at 1.25 × the floor (0.026 → 0.032 Å on dev-symmetry-lowering), with the
    same group in every case, and the shift tool names its lower group at that tolerance.
  - The split-site readings say `list_sites`, with κ₁, its z and δ ± its error.
  - The eight network cluster readings are gone (six records).
  - dev-jahn-teller's F sites 5 and 6 (κ₁ −0.95 and −0.93 at 216 copies, z −2.84 and −2.78) are no
    longer split. The F claim, from the site with the largest δ, stays.
  - dev-offcentre's claim reads δ 0.202 Å, where it read 0.201 (planted 0.2).
- The judgment subset's five A+ records move in the same ways, with every score unchanged. Their
  readings are in `arm-a-dev.json` and the `arms.test.js` pins, not here, since the held-out
  author may read this document.
- The moved A+ records are re-recorded in `arm-a-dev.json`; the A and D records stay as recorded.
- The test set was not run; its A+ records are re-recorded at the freeze.

**The protocol does not move.** What the model sees is byte-identical (the model-visible snapshots
pass without `-u`), so live runs under 2026-10-19.9 stay comparable. Arms A and A+ have no model,
and each of their records carries its own commit and protocol.

**The skills.** The skill lint gained a `rule` check (`RULES` in `skillLint/tables.js`). It flags
a skill sentence that states one of these rules otherwise. Five such sentences sit on its
allow-list with the skill's rule and the app's, awaiting the user's skill review (plan 5.3):
- displacement-anisotropy SKILL.md:27 and significance-reading.md:17: a split site is a negative
  Mardia kurtosis;
- coordination-shells SKILL.md:51: any polyhedron distortion goes to `analyze_clusters`;
- local-symmetry SKILL.md:38: a rung above the floor is a distortion;
- local-symmetry SKILL.md:59: whether a rung just above the floor counts is left to the user.

Approved by the manager on the user's delegation, 2026-10-19
(`validation/findings/decisions/2026-10-19-aplus-app-rules.json`; known answers in
`validation/__tests__/aplusRules.test.js`).

**Protocol 2026-10-19.10 (arm C without the method, plan item 1.8).** Arm C runs the turn with
`method: 'off'` ([Arm C](#arm-c-the-copilot-without-the-method)): no skills, no `EVIDENCE_RULES`, no
method sentence of `TOOL_RULES` or `ASK_RULES`, no method advice in a tool description or a check.
What arm C is sent changes: on `dev-offcentre` its system prompt is the `tools-ask.method-off`
variant (5,975 characters, the domain rules' inference parts gone too; the full method's `tools-ask`
is 8,575), its tool-free prompt is `no-tools.method-off`, its tool descriptions are
`harness/tools.arm-C.json`, a check that carried advice reaches it without
(`tool-results/*.method-off.json`), and the evidence gate's send-back names no analysis
(`harness/gate-send-back.arm-C.txt`). Arm B, the dock and arms A, A+ and D are unchanged: every
other snapshot file passes as it was. No arm C run exists under an earlier protocol. Each record
keeps the method its run followed (`method`).

**Protocol 2026-10-19.11 (the judgment subset's held-out matches, plan item 1.6 ii).** Five
judgment cases join the test set, one per dev judgment case, each with its situation, kind and
option: test-two-features (two_features: Jahn–Teller octahedra whose centre is also off-centred),
test-marginal-above and test-marginal-below (marginal: one weak zone-boundary order at one
amplitude, set from the spec's noncentral χ² so that the expected B power sits at 3σ after the
family-wise correction, two seeds), test-two-readings (two_readings: eight-site <111> off-centring)
and test-user-choice (user_choice: elongated octahedral clusters). A fresh author wrote them from
`validation/judgment-spec.md` alone: no dev judgment parameter, dev record, arm pin, finding,
transcript or PR description was opened, and no arm or model ran on them. Their values, the knob
explorations (on scratch seeds outside every case's sequence) and the seed rules' readings are
sealed with the other test cases since the holdout gate ([the sealed test
set](#the-sealed-test-set)); see `TEST_SET_SHA256`. Each case seed is the first of its own
sequence that a rule fixed beforehand accepts, from the case's engine readings, and each rule
checks the property the spec gives its situation:
- **two_features:** the centre split by `isSplitSite`, and the centre–ligand shell read as two
  peaks of CN 4 and 2. The first seed tried was taken.
- **marginal:** B's σ at q in the default `analyze_modulations` call (every element); the first
  seed in [3, 4] and the first in [2.5, 3), as for the dev pair. Both at rank 1 among the modes.
- **two_readings:** the centre not split and its PC1 unresolved, a clearly negative Mardia
  kurtosis, at the default tiling the map test at 3σ or more and the peak test below 3σ, the peak
  at 3σ or more near a <111> on a finer tiling, the centre–ligand shell 3 + 3. The first seed
  tried was taken.
- **user_choice:** the spec's four constraints (a deep trough, a joined shell at CN 4, the average
  edge beyond the trough, so that the first peak alone is a cut by the Workbench's CN margin, and
  an empty gap before the nearest contact between clusters). The first seed tried was taken.
- **The coincidence check.** The lead compared these values with the dev cases' without showing
  them, and two coincided: an A-site element (changed for another of the same spectator role and
  similar size) and an elongation (made smaller, with the thermal width lowered to keep the
  elongation-to-thermal ratio and the deep trough). Both rules were re-run unchanged and took their
  first seeds again; every reading above is the new configurations'.

`TEST_SET_SHA256` changes: results on the earlier test set are not comparable with this one.

**Protocol 2026-10-19.12 (the polarisation of an ordered mode, plan item 4.3; findings F2 and
F3).** What `analyze_modulations` returns changes for every arm that calls it:
- **The polarisation (4.3).** Both modulation engines give every reported point its
  `polarisation`: the leading axis of the real part of the element's cross-spectral tensor at q,
  as a Cartesian direction and a direct-lattice [u v w], with its share of the power
  ([ai-assistant.md, Step 19](ai-assistant.md#step-19--displacement-modulations-order-that-enlarges-the-cell-analyze_modulations)).
  The tool shows it for each ordered point (a `polarisation` row and its check: "polarised along
  [0 0 1] (5.5° off; 100 % of the power at q)").
- **F2.** A ratio or σ below its bar (3, 3σ) is never printed rounded onto it: a point at 2.98σ
  (the synthetic known answer) read "3σ" beside "≥ 3 is order" and a "no order" verdict; it now
  reads 2.98 and "below the 3σ bar".
- **F3.** The Šidák family is every element of the configuration, whichever the call names:
  `analyze_modulations({element})` corrected over that element's points alone, so naming one
  element of three took a point just below the bar over it. The note says so; a call without an
  element is unchanged.
- **Arm A+** reads a modulation claim's direction from the engine's polarisation at the q that
  decided the order, no longer from the PC1 of the element's most prolate site (a rule of its
  own, above), and records each ordered point's polarisation in its reading. A+ calls the engine
  for every element, so F3 changes none of its σ.

**Arms A and A+ on dev** (`arm-a --split dev` before at d9cb984, after on the code of 9a72a6f,
recorded under this protocol): every arm A and D record is identical, and every score. Two A+ records move,
each in its claim's direction and its reading, with the same score:
- dev-zone-boundary: the Ti claim reads [−0.05 0.09 1], 5.9° from the planted [0 0 1] as scored
  (the polarisation is 5.5°; the PC1 it read before, [−0.03 0.1 1], scored 6.0°: one site over
  512 cells, 1.5° of sampling error per component);
- one record of the judgment subset: its claim gains a direction, which its truth does not score
  (its values are in `arm-a-dev.json` and the `arms.test.js` pins, not here).

The two are re-recorded in `arm-a-dev.json`; the other records stay as recorded. The test set was
not run. What the model sees changes in the modulation tool's result only: its note
(`tool-results/analyze_modulations.json`) and a new sample of a call naming one element
(`analyze_modulations.element.json`); no tool description, prompt or first request changes, so the
32k budget pins hold. Approved by the manager on the user's delegation, 2026-10-19
(`validation/findings/decisions/2026-10-19-mode-polarisation.json`).

**Protocol 2026-10-19.13 (the tool transcript's stricter readings, architecture candidate #11).**
Every call now ends with one result status (`ok`, `refused`, `error`, `unknown_tool`, `skipped`),
decided by the tool contract, and every record step keeps it with the note the model was told and
the user's reply to a question ([Records](#records-transcripts-findings)). Two counts read them,
and what they mean changed:
- **"asked before"** in `unsafe_actions` counts only an earlier question that reached the user (an
  `ask_user` step with the user's reply as its `choice`), never a malformed or refused one, nor a
  window the Workbench settled with no card (it records no question);
- **number provenance** credits a number only to what the model was told: each tool message it
  read, its error and note included, the run context as shown and the question; a Workbench
  question's card, which the model never read, is no source.

What the model is sent is unchanged (every model-visible snapshot), and so are the safety counts'
other fields (`failed_calls`, `unknown_tool_calls`, now read from the status, and
`questions_asked`). Arms A, A+ and D carry neither count: on the dev set their records are
unchanged and are not re-recorded. The test set was not run. A live run under an earlier protocol
is not comparable with a later one on these two counts. Approved by the manager on the user's
delegation, 2026-10-19 (`validation/findings/decisions/2026-10-19-asked-before.json`,
`2026-10-19-provenance-told.json`).

**Protocol 2026-10-19.15 (a joined shell's split by number, plan item 4.2).** The review saw the
model subtract two g(r) peak positions to quote a Jahn–Teller split, a number provenance flags as
unsourced and no reference can carry. The shell catalogue now reads, for a shell joined across a
deep trough, its split: the parts cut at each inter-peak minimum (every part of a 2 + 2 + 2
octahedron), each with its CN (summing to the shell's), share and r²g mean distance, and long −
short with the configuration's sampling error
([ai-assistant.md, Step 9](ai-assistant.md#step-9--assembling-contextpair_correlations)). What
`find_gr_peaks` returns and says changes for every arm:
- **`first_shell.split`** (`parts`, `long_minus_short_A`, `long_minus_short_se_A`), only for a
  joined shell; the dataShape declares it.
- **The joined-shell note** quotes the parts' counts at their mean distances, long − short ± its
  sampling error and the shares by path, and no longer the peaks' separation: dev-jahn-teller's
  Cu–F reads 3.98 F per Cu at 1.929 Å + 2.01 at 2.251 Å, long − short 0.321 ± 0.0043 Å (planted
  0.32; the maxima of g are 0.31 apart).
- **The description** names the split: `dock/tools.json`, `harness/tools.arm-B.json` and
  `tools.arm-C.json` grow by 123 characters each (the budget pins: +48 tokens on the dock's worst
  turn, +49 on arm B's). The Demo joins no shell, so its `find_gr_peaks` results and every other
  model-visible snapshot are unchanged. The scripted runs' record files
  (`validation/__tests__/__snapshots__/records/`) move only in each request's length and hash
  (+121 characters a request with tools), and `harness.coordination-square` in its
  `find_gr_peaks` step on dev-jahn-teller's Cu–F shell (the split and the new note).
- **Arm A+** takes a Jahn–Teller claim's magnitude from the catalogue's split where it joins the
  claimed shell, and records the split as a `find_gr_peaks` reading; for a split shell the app
  does not join it keeps the separation of the maxima arm A reads, since the app gives no split
  there. Which shells are claimed stays A+'s own joining rule, and the long axis has no app rule
  (the g(r) is orientation-averaged), so A+ still names none. Arm A is unchanged.

**Arms A and A+ on dev** (`arm-a --split dev` before at c4d8d83, after at e9e67e1): one of the 45
records moves. dev-jahn-teller's A+ claim reads 0.321 Å where it read 0.31, with the split as a
new reading; its score stays 0.9167 (the magnitude was and is within its tolerance; the long-axis
family is still unnamed). Every arm A and D record, every other A+ record and every score is
identical. The record is re-recorded in `arm-a-dev.json`; the others stay as recorded. The test
set was not run. Approved by the manager on the user's delegation, 2026-10-19
(`validation/findings/decisions/2026-10-19-jt-split.json`).
**Protocol 2026-10-19.16 (the independent-cell null, finding F1; diffuse order between cells, plan
item 4.4).** What `analyze_modulations` returns changes for every arm that calls it:
- **F1.** Both modulation engines test an element's power against its own cells placed at random,
  which keeps every correlation inside a cell (the ligands of a Jahn–Teller octahedron, the
  vertices of a cluster) and destroys every correlation between cells
  ([ai-assistant.md, Step 19](ai-assistant.md#step-19--displacement-modulations-order-that-enlarges-the-cell-analyze_modulations)).
  The null of independent sites it replaces found order in 9–32 % of the configurations of the dev
  Jahn–Teller, two-feature and octahedral generator kinds (1,000 fresh ones each); the calibrated
  null finds it at the rate its bar allows. It also stops holding an independent site to the
  one-component bound: an isotropic site has three, so a weak order on one is read as far more
  significant than before.
- **4.4.** Each element gains `verdict.diffuseOrder`: a rank test of whether the power over the
  compared modes is spread in the pattern that correlation between neighbouring cells makes. The
  tool shows a `diffuse` row per element and says antiphase correlation as short-range order and
  in-phase correlation as neighbouring cells moving together (the Demo's three elements: in phase).
- **The marginal situation.** Both dev weak orders are now far above the 3σ bar, so the dev pair
  no longer straddles it; their amplitude was set from the null before F1. Re-setting it (and its
  held-out match's, set by the same rule) is the user's call; the case definitions are unchanged.
  (The dev pair was re-set at protocol 2026-10-19.19, below.)
- **Arm A+** reads cellDoubling, generalOrder and longWave as before; the diffuse finding is not
  among its readings (the verdict has no feature kind for it, and adding one changes the scoring).

**Arms A and A+ on dev** (`arm-a --split dev` before at c4d8d83, after on this protocol's code):
every arm A and D record is identical, and every score but one. One A+ record of the judgment
subset moves: its weak order is now ordered, A+ claims it and its score rises (its values are in
`arm-a-dev.json` and the `arms.test.js` pins, not here). It is re-recorded in `arm-a-dev.json`; the
other records stay as recorded. The test set was not run. What the model sees changes in the
modulation tool's result only (`tool-results/analyze_modulations.json` and
`analyze_modulations.element.json`: the σ of the zone-boundary rows, the note, the diffuse rows and
the Demo's in-phase note); no tool description, prompt or first request changes. The scripted runs'
record files that call it move with it (`fitCase.in-the-harness`, `harness.gate-split-site`,
`harness.over-claim`, `judgment.edge-above`, `judgment.edge-below`: the tool's rows and checks, and each
request's length and hash); the edge-of-significance harness test now finds its below-the-bar case at
a smaller amplitude by the weak pair's own rule, since both dev weak orders are above the bar. Approved by the
manager on the user's delegation, 2026-10-19 (`validation/findings/decisions/2026-10-19-modulation-null.json`,
`2026-10-19-diffuse-order.json`).

**Protocol 2026-10-19.18 (one cut rule: `compute_bond_angles` refuses on the CN margin).** Protocol
2026-10-19.8 left two readings of a window that cuts a shell, the Workbench asking on the CN margin
and the tool refusing on the coordination counts, with a recommendation for the user (below that
entry). The user approved it: one rule now decides, the shell module's `windowCut`
([ai-assistant.md, Step 9](ai-assistant.md#step-9--assembling-contextpair_correlations)). A window
that meets its first shell is cut when it holds fewer partner atoms than the whole shell by at least
max(0.5, 15 % of its CN), whatever its counts; the counts (more than half unbonded, under 60 % with
the modal count) decide only where the run reads no shell to compare: no partials file, no closed
shell for the pair, no density, or a window apart from the first shell (the module reads first
shells only). The Workbench asks on it and the tool refuses on it, unless `window_confirmed`. What
changes for every arm:
- **What the tool refuses.** `run.js shell-cuts` (now the Workbench's question and the tool's own
  refusal on every natural window, beside the counts rule): on the record's 149 windows (the ten
  dev cases it read and the Demo) the two agree on all 149 and 15 change verdict, as the readings
  predicted — 11 no longer refused (dev-offcentre's whole Nb–Nb shell, and ten windows with the top
  0.2 Å early, each within 15 % of its shell's CN), 4 refused now (the first peaks of the Demo's
  Ta–Se and Se–Se, dev-jahn-teller's four short Cu–F and dev-offcentre's one short Nb–O). With the
  judgment subset the dev and Demo windows number 227: 16 refused, the Workbench asking on the same
  16, 22 changed (16 no longer refused, 6 refused now).
- **What the model reads.** The refusal names the CN shortfall by path ("it holds 3.75 Se per Ta
  (shell_ab.cn_in_window) of the shell's 6 (shell_ab.cn, in shell_ab.window_A)…") instead of the
  share of central atoms; every result adds the shell reading (`shell_ab`, and `shell_bc` for a B–C
  window of its own) and the counts' `modal`, `modal_share` and `unbonded_share`, as information:
  no note says a count spread, since under a passing margin it was a whole shell running into the
  next or a window within 15 % of its shell's CN. The model-visible snapshots move in the five
  `compute_bond_angles` tool results (+197 characters of content for a whole shell, +248 for the
  refused spot-check window) and four Workbench tool messages; the tool's description and schema,
  the tool lists and the budget pins are unchanged. The scripted runs' record files move only in
  `harness.unsafe-window` (its step's data and the two later requests' lengths and hashes).
- **`facts.coordination.whole_shell`** of `compute_bond_angles` is the window's CN within the CN
  margin of the shell's, so a window apart from the first shell, or one holding more than it by
  the margin, is no whole-shell count for the coordination check.
- **The skills** still state the counts rule (coordination-shells SKILL.md lines 61, 78 and 79):
  the skill lint binds those numbers to the CN margin and holds them on its allow-list with both
  sides until the skill fact-sync PR.

**Arms A, A+ and D on dev** (`arm-a --split dev`, before and after): none of the 45 records moves
(A+ reads no bond-angle window); `arm-a-dev.json` stays as recorded. The test set was not run.
Approved by the user, 2026-10-19 (`validation/findings/decisions/2026-10-19-shell-cut-rule.json`).

**Protocol 2026-10-19.19 (the marginal situation re-set for the independent-cell null).** The
calibrated null of protocol 2026-10-19.16 gives an isotropic site three components where the
earlier null took one. The dev weak-order pair's amplitude had been set from that one-component
bound to an expected 3σ, so against the calibrated null both of its seeds sat far above the bar:
the situation was no longer marginal.
- **The knob** ([`judgment-spec.md`](../../web_app/frontend/validation/judgment-spec.md), the
  held-out author's only description of it). The B power at q is a scaled noncentral χ² of three
  degrees of freedom. The null is three weighted χ²₁ terms whose weights come from the site's own
  cell covariance. σ is Šidák-corrected over every zone-boundary point of every element. The
  expected significance is the median realised σ, which is the σ of the median power. The spec says
  how to evaluate it and how to check it with the engine on fresh configurations, and its scan
  (`judgment.test.js`) still finds no dev value. `judgment.test.js` evaluates the formula with the
  engine's own functions (the noncentral χ² median against scipy's).
- **The dev pair.** The pair's amplitude is re-set by that formula to an expected 3σ. Over 5,000
  fresh configurations the engine clears the bar on about half. The seeds are re-taken by the
  pair's rule, fixed before any σ was read at the new amplitude: the same sequence, the first
  realised in [3, 4] and the first in [2.5, 3). Every other parameter is kept. The values are in
  the decision record and the PR, not here, since the held-out author may read this document.
- **Arm C.** `analyze_modulations`' long-wave and in-phase notes lose their causal parentheticals
  ([Arm C](#arm-c-the-copilot-without-the-method)). Arm B and the dock read them byte for byte as
  before. Arm C's samples are `tool-results/analyze_modulations.method-off.json` and
  `analyze_modulations.element.method-off.json`. On the dev set three cases carry a long-wave note
  and none an in-phase note.

**Arms A, A+ and D on dev** (`arm-a --split dev` before at f0979a7, after at c9fea9c): the six
records of the two re-set cases move, with their configurations and truths. One A+ record claims
the weak order as before. The other no longer does, since its order is below the bar, so A+ misses
there as the situation intends; its score falls. Every other record is identical apart from
protocol and commit, and so is every arm A and D score. The six are re-recorded in
`arm-a-dev.json`. The test set was not run; its two marginal matches are to be re-set from the spec
by their author, in a separate change. What the model is sent is unchanged for arm B and the dock
(every model-visible snapshot). The scripted record `judgment.edge-above` moves with its case (the
tool's rows and check, the truth, and each request's length and hash); `judgment.edge-below` keeps
its configuration and does not change. Approved by the manager on the user's delegation,
2026-10-19; the user may reverse it before the freeze
(`validation/findings/decisions/2026-10-19-marginal-reset-dev.json`).

**Protocol 2026-10-19.20 (the held-out marginal pair re-set).** The marginal situation was re-set
after F1's corrected null raised the test's power, 2026-10-19. A fresh author, who had seen no dev
value, re-chose the two held-out marginal matches' amplitude afresh from the updated
[`judgment-spec.md`](../../web_app/frontend/validation/judgment-spec.md), with new seeds; every other
parameter is kept. The values are sealed with the other test cases (see `TEST_SET_SHA256`).
- **What the author was given:** the spec, the generator code (`validation/cases.js`), the sealed
  module it edits and its loader, and the modulation engines to measure its cases. Nothing else:
  no dev judgment parameter or record, no arm pin, no decision record, no doc, no git history.
- **How:** the amplitude from the spec's formula alone, checked with the engine on many fresh
  configurations (about half clear the bar), and each seed by a rule written down before any σ was
  read: the first of a fixed sequence realised in [3, 4] for the match above the bar, in [2.5, 3)
  for the one below.
- **Blindness, checked from the author's tool log by the lead:** every file it opened is in the list
  above, and no dev value reached it. One disclosure: the sealed module's comment on protocol
  2026-10-19.11 named the two test values that had coincided with dev values; the author read it.
  Those concern two other judgment cases, not this pair, and the comment no longer names them.
- **The holdout scan** no longer counts the digits that end an identifier as a number: the new
  seed matched inside `TEST_SET_SHA256` beside the word "seed" on this page.
- **Arms A and A+ on test** are re-recorded with `--final` from a clean tree in the sealed folder.
  Approved by the manager on the user's delegation, 2026-10-19; the user may reverse it before the
  freeze. `TEST_SET_SHA256` changes: results on the earlier test set are not comparable with this
  one.

**Protocol 2026-10-24.1 (the fit-assessment interview's code half).** The user's answers in the
rmc-fit-assessment interview, 2026-10-24, change what the model is sent on a case with fit files
(`dev-fit-unconverged`; the Demo) and what the checks on the answer read, never a case, the question,
the scoring or the verdict schema:
- **Convergence on every χ² term and the total.** `assess_convergence`'s description and result,
  the context's `convergence` section (`last_10pct` per term and for the total, in place of
  `change_over_last_10pct` and `recent_slope_per_step`) and the system prompt's `watchdog_status`
  sentence read one verdict: converged only when no χ² term of the .log and not their total still
  falls beyond 3 × its noise over the last tenth (ai-assistant.md, Step 14). On the fit case the
  status stays `improving` (both terms and the total fall); its truth is unchanged. The context
  budget is 6,500 characters (6,000 before).
- **The constraint-limited rule:** `find_gr_peaks` calls a peak constraint-limited above 5 % of the
  peak at the cut; no dev case has a constraint, so no dev record reads it.
- **Two checks on the answer,** `long-wave` (order claimed at a box-scale wavevector) and
  `single-configuration` (a distribution's numbers quoted without saying they come from one
  configuration): any finding lowers a run's computed Outcome to "not verified", so arms B and C
  records under an earlier protocol are not comparable on the Outcome. `audit_eval` still maps a
  spurious modulation claim to the `superstructure` rule only, so a `long-wave` finding on one shows
  among `unmatched_findings`.
- **Arms A, A+ and D on dev, before (8edf8c0) and after:** all 45 records identical; the baseline is
  not re-recorded. Approved by the user, 2026-10-24
  (`validation/findings/decisions/2026-10-24-fit-assessment-interview.json`).

**Protocol 2026-10-24.2 (the interview answer checks of four skills).** The user's answers in the
coordination-shells, local-symmetry, displacement-anisotropy and local-order interviews, 2026-10-24,
change what the model is sent and what the checks on the answer and the evidence gate read, never a
case, the question, the scoring or the verdict schema:
- **What the model reads.** `EVIDENCE_RULES` names `polyhedral_distortion` for distorted polyhedra;
  `analyze_clusters`' method keeps it to finite clusters (arm B; arm C reads no method);
  `find_gr_peaks`' description and result give the average structure's first shell for the pair
  (`first_shell.average_structure_d_A`, with a warning past the first peak's FWHM);
  `analyze_displacement_directions`' description and result label the peak a finding, a hint or
  nothing (`peak.tier`, the `neighbour` block for a candidate); `analyze_modulations` says in-phase
  correlation as "neighbouring cells moving together (correlated motion), not order". The symmetry
  tools' results are unchanged.
- **Checks on the answer:** `window`, `subgroup`, `cell-doubling`, `disorder`, `no-order`,
  `correlated-motion` and `no-correlation` are new; `coordination` reads "n A per B" as B's count
  and a clause stating the whole shell's window as a claim, and `direction` needs a finding of the
  tiers (ai-assistant.md, *Checks on the answer*). Any finding lowers a run's computed Outcome to
  "not verified", so arms B and C records under an earlier protocol are not comparable on the
  Outcome. `audit_eval` still maps a spurious modulation claim to the `superstructure` rule only, so
  a finding of the new modulation checks on one shows among `unmatched_findings`.
- **The evidence gate:** distorted polyhedra are met by `polyhedral_distortion` where the turn offers
  it, else `analyze_clusters` (an `anyOf` row asks for the first analysis the turn offers); a split
  or flat-topped site called disordered goes back for `analyze_modulations`; "no order" also needs
  `polyhedral_distortion` on each centre of polyhedra no finite cluster holds, where offered (no dev
  run offers it yet).
- **Arm A+** follows the tool's direction tiers and reads the diffuse verdict (antiphase: a
  modulation claim with no q; in phase: read, never claimed). **Arms A, A+ and D on dev, before
  (afa45e1, the branch's base) and after each of the four changes:** no verdict or score moves. Four A+
  records' readings add the tier (dev-offcentre, dev-jahn-teller, dev-jt-off-centre,
  dev-octahedral-clusters: every split site's peak is a finding), and every A+ record's
  `analyze_modulations` reading adds the diffuse verdict (no dev case passes it). The baseline is not
  re-recorded: its verdicts and scores stand. Approved by the user, 2026-10-24
  (`validation/findings/decisions/2026-10-24-coordination-shells-interview.json`,
  `2026-10-24-local-symmetry-interview.json`, `2026-10-24-displacement-anisotropy-interview.json`,
  `2026-10-24-local-order-interview.json`).

**Protocol 2026-10-24.3 (a new tool, `polyhedral_distortion`).** The user's interview answer 7 for
the coordination-shells skill: corner-sharing or extended networks of polyhedra are judged by angle
spreads plus a standard polyhedral distortion index, never by the cluster null. Arms B and C are now
offered `polyhedral_distortion` (Baur's DI, Robinson's quadratic elongation and bond-angle variance
of a centre's coordination polyhedra, over the configuration and the average structure;
[ai-assistant.md](ai-assistant.md), Step 21). It refuses a window that cuts the shell by
`compute_bond_angles`' rule, and arm C reads its description without the method advice. Arms A and
A+ do not read it: their 45 dev records are identical apart from protocol and commit, and the
baseline is not re-recorded. The test set was not run. What the model is sent changes by the tool's
definition (arm B +1,741 characters, arm C +1,246) and the first request's tool list. The scripted
records add the tool to `tools_offered`, and each request's length and hash move. Approved by the
manager on the user's delegation, 2026-10-24
(`validation/findings/decisions/2026-10-24-polyhedral-distortion.json`).

**The freeze (2026-10-24).** The validation is frozen at commit 226f6f2 (protocol 2026-10-24.4;
`TEST_SET_SHA256` as pinned in `split.test.js`), recorded in `validation/baseline/freeze.json`. Arms
A, A+ and D were re-recorded with `--final` from a clean tree at that commit (`arm-a --split dev
--final`, `arm-a --split test --final`): on dev A 0.814, A+ 0.952 and D 1.000 over 15 cases; on test
A 0.800, A+ 0.948 and D 1.000 over 14. No score differs from the baselines before it; the records
gain only the readings later changes added (the direction tier, the diffuse reading). No tag: tags
are the user's. The live arms B and C wait for the user's go/no-go.

**Protocol 2026-10-24.5 (a refused analysis backs no claim, the user's decision A5).** The user's
words: "A5: option 2. A refused analysis doesn't count: the gate sends the claim back, and the model
re-runs on the whole shell or asks me to confirm the window. Land it before the B/C runs." The
evidence gate counts only a call that answered: a `refused` one (the tool's own `fail` check, such as
`compute_bond_angles` refusing a window that cuts the shell) backs no claim, under `anyOf`,
`firstOffered` and `allOf` alike; until now it ran and backed its claim. It is still an analysis (the
dock's footer, each record step's `analysis`, `refused_calls`) and the Outcome still reads it as a
failure. The send-back for a claim a refusal would have backed says it refused: arm B's names the
tool and, for a window cut (a refusal whose data carry `refused`), both fixes — the whole shell, or
the user's own window through `ask_user` and `window_confirmed` —, any other refusal only "Run it so
that it answers"; arm C's names no tool, window or remedy. The end-of-turn finding and the Outcome's
reason name the refusal ("did not run or refused"). The system prompt, the tools and every
model-visible snapshot are unchanged (the 48k headroom stays 4,241 tokens for the dock's worst turn
and 7,044 for arm B's). Arms A, A+ and D run no model and no gate: on dev their 45 records are
identical apart from protocol and commit. The test set was not run; one re-freeze follows the next
change. Approved by the user, 2026-10-24
(`validation/findings/decisions/2026-10-24-refused-not-evidence.json`).

### The report beside the case score (plan item 1.6, part i)

The engines alone (A+) score 0.99, because most components of a case are facts they get right every
time: a perovskite case has four symmetry fields and the spurious-claim test beside one planted
feature. On the case score alone B can at best tie A+ ([the bar after batch 1](#arm-a-results-no-model-2026-10-14-commit-b4fe09f)).
So under each split's tables, `report` takes the case score apart (`comparison.js`, from the run
records). It adds three sections:

1. **The case score's parts.** For each column, the case score restricted to one group of its
   components (`scoring.js` `scoreComponents`, `restrictedScore`):
   - the **feature-only score**: the planted-feature components alone (1 found and accurate,
     0.5 a quantity off, 0 missed). A case with no planted feature has none and is left out, so
     the cell says over how many cases it is taken;
   - **no spurious claim**: the share of runs with no claim beyond the planted and tolerated ones;
   - **the symmetry fields** (a case whose truth defines none is left out).

   Each case weighs once. Under each live column, a row gives its difference from A+, taken on
   the cases both have.
2. **Component by component, against A+.** Every component of every case, with arm A's value,
   A+'s, each live column's (its mean over the runs) and the live column minus A+. A non-zero
   difference is bold, and a line per live column counts the components higher, lower and the
   same.
3. **Every disagreement with the planted truth.** For each case, each feature or fact a run got
   wrong, with the answer, the truth (with its tolerance) and the columns that made it, with
   their runs (2/3):
   - a missed feature, with what took its place: another kind on its element ("wrong kind",
     an ordered mode read as off-centring) or its kind on another element ("wrong element");
   - a direction, magnitude or q off;
   - a spurious claim;
   - a symmetry field wrong or unanswered;
   - no parseable verdict.

   Tolerated claims (true partial descriptions) are not disagreements. Each disagreement costs a
   component, and each component below 1 has at least one. The cases where no run disagrees are
   listed after the table.

**How to read it.**
- B − A+ on the feature-only score is what the model adds or loses on the planted order.
- A lower "no spurious claim" is an over-claim the engines did not make.
- A disagreement that A+ and B share is the engines' limit (the JT long axis on dev). One that
  B makes alone is the model's own. One that A+ makes alone is where the model added something.

**On the dev baselines** (`arm-a-dev.json`; no B run yet, so the report says so and shows A
beside A+). These are the ten dev cases outside the judgment subset; the subset's results are kept
out of this document ([the judgment subset](#the-judgment-subset-plan-item-16-part-ii)), and
`comparison.test.js` pins this table on the ten:

| Column | Case score | Feature-only score (6 cases) | No spurious claim | Symmetry fields (9 cases) |
| --- | --- | --- | --- | --- |
| A | 0.825 | 0.25 | 0.90 | 0.92 |
| A+ | 0.99 | 0.92 | 1 | 1 |
| D | 1 | 1 | 1 | 1 |

A+ has one disagreement: the Cu long axis in dev-jahn-teller is not given (A has the same one).
A's other disagreements:
- dev-zone-boundary: modulation Ti missed, read as off_centre Ti (wrong kind), which is also a
  spurious claim;
- dev-breathing: cluster_breathing Mo missed;
- dev-cell-wrap: rigid_cluster V missed;
- dev-symmetry-lowering: symmetry_lowering Ta missed; Pm-3m given for P4mm, and called
  centrosymmetric and non-polar.

On dev, B can gain over A+ at most 0.08 on the feature-only score (the JT axis is half of one of
six features), against 0.008 on the case score.

**What does not change.** The case score, every existing export of `scoring.js` and the arm A/A+
pins (`arms.test.js`) are unchanged. `scoreComponents` reads the components back from
`scoreVerdict`'s result, and `__tests__/comparison.test.js` checks three things on every dev
baseline record, the truth and no verdict: their mean is the case score; the baseline re-scores
exactly as recorded; and the listing matches the components. **Not a protocol change.** The report
reads records: it changes no record, no score and nothing a run sends, so `PROTOCOL_VERSION`
stays 2026-10-19.7, and records under every protocol go through it. Plan item 1.6's part (ii),
the judgment cases, is a separate change (protocol 2026-10-19.9).

## Records, transcripts, findings

Each live run writes:
- `<case>__<arm>__<model>__r<n>.json`: the view the run was opened with and what the model was
  told about it (below), every step with who decided it (the model or the Workbench), the
  reasoning and text the model wrote before the call, the arguments, the checks and the result it
  read, then the final answer, verdict, score, answer checks, provenance, safety, and the truth;
- `.md`: the same as a self-contained transcript, with numbered steps and the planted truth last.

**What the run had (plan item 1.3).** The app's Copilot can be handed more than a case gives the
harness, so each record states it (and, since protocol 2026-10-19.10, `method`: the method the turn
followed, `'full'` in arm B and `'off'` in arm C):
- `view`: the run module's normalised view, `{ toleranceA, hiddenPlots, pcaSites, liveData,
  ladder }`. By default it is `HARNESS_VIEW`: 0.2 Å with the ladder, nothing hidden, no PCA table,
  Live Data off. In the app the tolerance is the user's ladder pick, and the PCA table joins once
  the PCA Ellipsoid page has computed it;
- `context`: `sections`, the context sections the model read (`harness.js` `CONTEXT_SECTIONS`, in
  `buildRunContext`'s order); `absent`, the sections it lacked; `earlier_turns: 0`, since every run
  is one fresh question (the app sends up to eight earlier turns). A case has `structure`,
  `symmetry`, `pair_correlations` and `datasets` (its partial g(r)). It has no
  `configuration_optimization` (no move counters), `run_settings`, `partial_weights` or
  `convergence`, and no `pca_displacements` unless the PCA table is on. The exception is
  dev-fit-unconverged, whose fit files add `configuration_optimization`, both fitted datasets
  with their Rw, `partial_weights` and `convergence` (it has no `run_settings` either). Read from
  the context as the model was shown it (the turn's `start`);
- `requests` (since protocol 2026-10-19.4): every model request as it was sent, in order, as
  `{ round, why, sent }`, and since 2026-10-19.5 `num_ctx`, the context the guard sent Ollama for
  it (null on other servers; since 2026-10-19.7 read from `sent.plan.numCtx`): `why` is `first`, `tools` (the
  results of the round before), `reminder`
  (after an empty or looping reply), `gate` (an answer the evidence gate sent back), `final` (the
  analysis budget used: `tool_choice` 'none') or `no_tools` (the retry of a model that refuses
  tools); `sent` is the request as the port prepared it (`providerModel`: base URL, model id,
  temperature (no key at all for a model that takes none, since 2026-10-24.6), the reply cap
  `maxTokens` where the model sets one (since 2026-10-24.6), messages, tools, tool choice and,
  since 2026-10-19.7, `plan`, the context guard's
  decision: `server`, and on Ollama `numCtx`, `numPredict`, `estimate`, `needed`,
  `loadedContext`, `maxContext`, on LM Studio `estimate`, `needed`, `loadedContext`, with the tool
  choice's form `toolChoice`: `body`, `omitted` or `tools_withheld`; never the API key). On an
  Ollama final request `tools` is absent (withheld). A request the spend cap or the context guard
  refused was never sent and is not there; the run's `error` names the refusal (since
  2026-10-19.7 for the guard, which decides in `prepare`). A request sent and then stopped (Ollama
  read only part of the prompt; a model LM Studio loaded for the request at too small a window)
  is there, followed by the error.
- `stopped_by` (architecture candidate #10; no protocol change, since nothing sent moves): how
  the turn ended, from its terminal event (`done`, `stopped` or `failed`), as one of the turn
  policy's `STOPPED_BY` keys (`src/llm/copilot/turnPolicy.js`): `answer`, `max_rounds` (the
  analysis budget or the turn's time used), `no_progress`, `empty`, `stopped` (Stop) or `error`
  (a stall included). The dock's note under a turn and the computed Outcome's "cut short" read
  the same table. A record from before it has no `stopped_by` key; its `final.stopped_by` (an
  answered run) is the same value.

**Each step (architecture candidate #11, the tool transcript).** A step is the call's `tool_end`
as the loop sent it, with who decided it: `step`, `round`, `id`, `name`, `by` (`model` or
`workbench`), `args`, `why`, `status` (the severity: `ok`, `warn` or `failed`), `result_status`
(how the call ended, the tool contract's: `ok`, `refused`, `error`, `unknown_tool` or `skipped`),
`checks`, `data` (as the model read it; for a Workbench question, its card's reply, which the
model never read), `error`, `note` (what the app added for the model beside the result; for a
Workbench question, the text the model read in its place), `choice` (a question the user
answered: `source`, `index`, and `confirmed`, whether a window question's answer keeps the user's
own window), `ms`, `cached`, `ref` (null with an error and for a Workbench question) and
`analysis`. A step still running when the turn stopped or failed is listed too, `status:
'stopped'` with `result_status` null, and no count reads it. `safety` gains `refused_calls`. These
keys were added without a protocol change (nothing sent moved); the readings that use them came
with protocol 2026-10-19.13 (*Reported beside the score*).

**What the run ran under (plan item 1.2).** `policy` (since protocol 2026-10-19.5): the run policy,
the server and the window the run ran under, the limits, the temperature, the reply cap, the context
guard's decisions and the final round's form ([The run policy](#the-run-policy-plan-item-12)).

The transcript prints the view, the context and the policy under its header. **The PCA-table option** (`run.js live --pca-table`,
dev set only; `runLiveCase({ pcaTable: true })`) gives each run the PCA table from its own engine,
the rows `list_sites` reads, as after a visit to the PCA Ellipsoid page. It is a dev check of what
that section changes and is off by default. Such runs write to their own folder (`…_pca-table`),
and `report` gives them their own column (`<model> · <arm> + PCA table`). It changes the context the
model reads, so those runs are not comparable with the default runs.

`report` puts every record side by side (case × model·arm), keeps failed runs as score 0, takes
each split's case score apart against A+ ([the report beside the case
score](#the-report-beside-the-case-score-plan-item-16-part-i)), sets the columns' run policies
side by side ([what a comparison across models also
compares](#what-a-comparison-across-models-also-compares)), and lists every run used.

Supervisor findings (`findings.js`, `run.js findings <file.jsonl>` validates them) take the form
`{case, run, step, kind: bug | missing-tool | skill-gap | orchestration | science-error, evidence,
proposed_fix, needs_user_approval}`:
- bugs are fixed with a regression test and merged on green CI;
- `missing-tool` and `skill-gap` always need the user's approval (batched per cycle);
- a test case may only report a bug.

Where a finding stands is recorded beside it: `status` (`proposed`, the default, `approved`,
`implemented`, `rejected`), `approval: {by: 'user', date, batch}`, the `commit` that implemented it
and the regression `tests` that pin it. An approved or implemented finding that needs approval
records it; an implemented one names its commit, and an implemented code change names its tests.
Batch 1 (`validation/findings/batch-1.jsonl`) holds every tool change of the
[review](copilot-review.md), approved by the user on 2026-10-14 and implemented on 2026-10-17:
L1–L5 (L4 is two analyses) and the five skill proposals, eleven findings (the skill ones came with
the skills PR, which waits for the user's review). `split.test.js` validates the file.
Decisions on the validation itself take the same form, one JSON file per record in
`validation/findings/decisions/`, named by the record's id (the date and a slug,
`2026-10-19-list-sites-split.json`), so that parallel changes never edit the same file.
`parseFindingFiles` (`findings.js`) reads every `*.json` there in id order and validates each by the
findings schema, in `split.test.js` and in `run.js findings <folder>`. So far four, all approved by
the manager on the user's delegation on 2026-10-19:
- the recomposed test breathing case (protocol 2026-10-19.2);
- plan item 4.1, `list_sites` showing a split site's static displacement
  ([ai-assistant.md](ai-assistant.md), Step 20);
- the engines' own verdicts (architecture candidate 5);
- the judgment dev subset (plan item 1.6, part ii; protocol 2026-10-19.9).

The schema's approver is the user, so the delegation is stated in each record's evidence and
proposed fix.

### The answer corpus (plan item 3.5)

The evidence gate and the checks on the answer (the verification,
`src/llm/copilot/verification/`: the claim and facts catalogues) read claims with regexes.
`src/__tests__/answerCorpus.test.js` runs them over real answers and records what they read:
- **The corpus:** one file per answer in `src/__tests__/answerCorpus/entries/`. Each file holds:
  - the text and its `fidelity`: `verbatim` (committed as the model's words), `transcribed` (read
    off a committed recording) or `quoted` (a docs or skill quotation, possibly shortened);
  - its `source` (file, lines, commit), the model, the question and the run (the Demo or a dev
    case, never a test case);
  - the steps that ran, by name and arguments, or `null` when the record does not show them.
- **The steps run** with the real engines. They must reproduce the numbers the transcript shows
  (`reproduces`).
- **One adapter** turns an answer into the pinned reading (`readAnswer` in
  `answerCorpus/readAnswer.js`, through the verification's interface only: the steps recorded,
  then the answer judged):
  - the claims the gate reads;
  - the gate's send-back given the steps, or "steps unknown";
  - the checks that fire.
- **Each reading is judged by a person.** Every send-back, every finding, or the silence carries
  `right` and a one-line `reason`. The test lists the false alarms and the misses.
- **When one verification module replaced the gate and the checks** (architecture candidate
  #4, 2026-10-19), only `readAnswer` changed (and its two call sites): the pins stayed.

**On 2026-10-19** the corpus holds nine answers on the Demo:
- **Verbatim:** gemma4:26b's three correct answers of 2026-10-11.
- **Transcribed:** Muse Glimmer's GaSe₄ answer from the README GIF.
- **Quoted:** five wrong answers from the 2026-10-10 spot check and review.

Four readings fire, and all four are right:
- chirality and inversion lost, on qwen3.8's "P-42₁m is chiral and centrosymmetry is lost";
- convergence, on "Converged — a plateau.";
- the gate's cluster send-back, on "the cluster is distorted" with no `analyze_clusters`.

There is no false alarm and one miss: "no cluster distortion", concluded from a pooled direction
map, needs no analysis under any rule.

No real answer yet makes a claim that the gate's split-site, lowering, modulation, direction or
no-order rules read, or that the superstructure, cluster, polar or coordination checks read. The
live runs' records are what will exercise them.

**Adding an answer:**
1. Copy the final answer of a live run's record on a dev case (or of a Demo session) verbatim into
   `entries/<date>-<model>-<topic>.json`. Never use a test case's record.
2. Copy the steps from the record, and a few of its numbers into `reproduces`.
3. Run the test. Its diff gives the reading to put in `expected`.
4. Judge each send-back and finding, or the silence.

**A change to a claim regex** must keep every pin, or update the entry in the same PR with a
judgement for each new reading. A new false alarm or miss shows in the test's lists.

## The run policy (plan item 1.2)

The cases, the question, the prompt, the tools and the scoring are the same for every model. The
rest of what a run runs under is its **run policy** (`validation/runPolicy.js`). Each record keeps
it (`record.policy`), each transcript prints it under its header, and the report sets the columns'
policies side by side.

| Field | What it states |
| --- | --- |
| `server` | The kind the context guard detects (`detectServer`: `lmstudio`, `ollama`, or `openai` for any other server, the cloud providers among them), Ollama's version, and the kind the model expects (`models.js` `server`) |
| `temperature` | The model's (0.2), and the floor of a request asked again after a looping reply (0.6). Claude Sonnet 5.5 takes none (since protocol 2026-10-24.6): `{ base: null, warmer_retry_min: null, sent: false }`, no temperature key in any of its requests, the retry included |
| `limits` | The turn's budgets as the run used them (`runPolicy.js` `TURN_LIMITS`: the agent loop's `DEFAULT_LIMITS` and the verification's send-backs, `maxEvidenceGates`): 6 rounds of at most 4 calls, 10 minutes, 3,000 characters a result and 12,000 a turn, 2 questions, 2 gate send-backs, the stall and no-progress guards |
| `context` | The minimum window for the server kind; how the window is sized; the window read before the run (and after it, on LM Studio); the model's maximum (Ollama); every `num_ctx` an Ollama request was sent (`num_ctx_used`; each request in `requests` keeps its own) |
| `reply_cap` | Ollama: `num_predict`, the num_ctx minus the prompt estimate. LM Studio: none sent (the guard needs 4,096 tokens of room). Cloud: none sent, unless the model carries a cap: Claude Sonnet 5.5 sends `max_tokens` 8192 in every request (since protocol 2026-10-24.6; `{ sent: 'max_tokens', max_tokens: 8192 }`) |
| `context_guard` | The guard's decisions: a request that would not fit is refused, never cut; 2.5 characters a token, 1.4 for tool JSON, 8 tokens a message; 8,192 tokens reserved for the reply; 4,096 of room on LM Studio; contexts in steps of 8,192, at least 16,384; the 0.4 read-share backstop |
| `final_round` | What the last request sends for `tool_choice: 'none'` ([below](#the-final-rounds-tool_choice)) |
| `history` | One fresh question; earlier reasoning is never re-sent |
| `spend` | An API model's cap ($20 for Claude Sonnet 5.5), and the reply each request is priced at (8,192 tokens, the cap it is sent with) |
| `warnings` | A window that changed during the run |

**The decision (2026-10-19, the manager on the user's delegation).** The context guard, its
estimate, the 8,192-token reply reserve, the loop's budgets and the prompt stay as they are:
refusing a request rather than cutting it is the behaviour wanted. Instead, the live runs give
each model a window that holds the worst turn (`MIN_CONTEXT_TOKENS`):

| Server | Minimum window | Sized by |
| --- | --- | --- |
| LM Studio | **49,152 tokens (48k)** | The window the model is loaded with; read before and after each run |
| Ollama | None | The guard, per request (`num_ctx`: at least 16,384, at most the model's maximum); recorded |
| Cloud | None | The provider; never probed |

48k holds the worst realistic turn of the dock and of arm B (two skills read whole, a full round of
analyses, the 8,192-token reply reserve; `src/__tests__/contextBudget.test.js`), with 10,452 and
13,182 tokens to spare (arm B's worst is dev-fit-unconverged since protocol 2026-10-19.6; protocol
2026-10-19.8 added 90 tokens to the dock's and 16 to arm B's, protocol 2026-10-19.15 48 and 49).

**The pre-run check.** `run.js live` calls `preflight` before every run. It reads the model's
server the way the guard does: `detectServer`, then LM Studio's `GET /api/v0/models`
(`loaded_context_length` of the entry whose id is the model's). It refuses the run when:
- the server at the model's URL is another kind, or nothing answers there;
- an LM Studio model is not loaded under its id (a request would load it at LM Studio's default);
- an LM Studio model is loaded below 48k.

Nothing is sent and no record is written, so the same command resumes once the model is
reloaded. The message says what was found, what is needed and how to reload:

> LM Studio runs qwen3.8-27b with a 32,768-token context, but the run policy needs at least 49,152
> (48k) for an LM Studio model, so the run was not started. In LM Studio, eject qwen3.8-27b and
> load it again with a Context Length of 49,152 or more (its load settings), or from a terminal:
> lms unload qwen3.8-27b && lms load qwen3.8-27b --context-length 49152. Then run the same
> command again: the runs already recorded are kept.

After the run the window is read again. If it changed (the model reloaded, or loaded again by a
request at LM Studio's default), the record, the transcript and the report carry a warning.

**Claude Sonnet 5.5's requests (protocol 2026-10-24.6; the user's decisions, 2026-10-24: "Fix
both first").** Sonnet 5.5 takes only its default sampling: a request carrying a temperature is
refused (HTTP 400), so no run could have started. Its `models.js` entry sets `temperature: null`
(as the sibling project NEBULA does) and `maxTokens: MAX_REPLY_TOKENS` (8,192). The model port
(`providerModel`) then sends no temperature key in any request, the warmer retry included (the
retry is still asked, unchanged), and `max_tokens: 8192` in every body, so a reply, hidden thinking
included, stops where the spend tracker's worst case already priced it. Both are settings, never
read from the model's name: the dock passes neither and its requests are as before, and the local
models' requests are byte for byte as before (`runPolicy.test.js`). `max_tokens` is never sent by
default (OpenAI's newer models refuse it).

**The dock at 32k.** The decision fixes the validation, not the dock. A user who loads
qwen3.8-27b in LM Studio at 32k can have the worst two-skill turn refused. The budget test pins the
dock's worst turn 1,836 tokens over by LM Studio's rule (5,932 by `planContext`'s), and arm B's
3,202 over by `planContext`'s. The fix is on-demand skill references (plan item 5.5, second half):
the skills' long reference material leaves what `read_skill` returns whole. That is skill text, so
it waits for the user.

### The final round's tool_choice

When the analysis budget is used, the loop asks for the answer with `tool_choice: 'none'`. The
model port puts it in each server's form (`provider/model.js` `toolChoiceForm`; the request as sent
names it, `sent.plan.toolChoice`). The scripted check (`__tests__/runPolicy.test.js`) sends a final
round through the production port (`providerModel` → `streamCopilotChat` → `client.js` or
`ollama.js`) to each provider preset, with a fake server, and pins what goes out:

| Provider | Endpoint | The final request sends | Form | Honoured? |
| --- | --- | --- | --- | --- |
| OpenAI | `/v1/chat/completions` | `"tool_choice": "none"`, the tools still sent | `body` | Documented by OpenAI |
| Gemini (OpenAI-compatible) | `/v1beta/openai/chat/completions` | The same | `body` | Not verified (no live call) |
| Anthropic (OpenAI-compatible) | `/v1/chat/completions` | The same | `body` | Not verified (no live call) |
| LM Studio | `/v1/chat/completions` | The same | `body` | Not verified; the smoke run shows it |
| Ollama | `/api/chat` (native) | **No tools** and no `tool_choice` (the native API has none) | `tools_withheld` | **Enforced** (since protocol 2026-10-19.7): with no tools the model cannot call one |

Before protocol 2026-10-19.7 Ollama's final request still offered the tools, and a reply of calls
only ended on the loop's empty-final answer. Every other round on Ollama still sends the tools and
no `tool_choice` (form `omitted`: its default, auto). The loop keeps its guard for the servers
that take the field: calls in a final reply are never run, and a final reply with no text ends on
the fixed "budget used" answer.

### What a comparison across models also compares

The report adds a table with one row per column: the server, the window, the reply cap, the final
round's form, the temperature, the limits and the policy warnings. It says whether the limits, the
temperature and the context guard are the same in every run, and lists what no policy sets. Today
a comparison across the four models also compares:
- **the context sizing:** Ollama gets a sized `num_ctx` and a `num_predict`; LM Studio a loaded
  window of 48k or more and no reply cap; the API model its own window and `max_tokens` 8192;
- **the temperature:** 0.2 for the local models (0.6 or more on the warmer retry); none sent to
  the API model, which takes only its default sampling;
- **the final round:** `tool_choice: 'none'` on the OpenAI-compatible servers, honoured where the
  server honours it, and no tools on Ollama, which is enforced;
- **what no policy sets:** sampling beyond the temperature (top_p, top_k, min_p, the repetition
  penalty), the reasoning mode and its effort, the tokenizer (the guard's estimate of a request is
  the same for every model, the real count is not), the quantization and the server's release.

## Models and spend

Each model gets the same cases, prompt, tools and scoring, and the results are reported side by
side, each under its run policy:
- qwen3.8-27b (LM Studio, loaded with a context of 48k or more);
- gemma4:26b (Ollama; the guard sizes `num_ctx`);
- meta/muse-glimmer (LM Studio, 28B, tool use, up to 128k; loaded with 48k or more);
- Claude Sonnet 5.5 (API).

The API model (`claude-sonnet-5-5`) runs the test set only, arm B, 3 runs per case, under a $20
cap, RMC's share of the cycle's $40 total, which may be up to $28 (the user's decision,
2026-10-24: "Raise to $20"; $12 before). The spend tracker prices each request before sending it (prompt
estimate + 8,192 tokens of reply, the `max_tokens` each request is sent with since protocol
2026-10-24.6, so that worst case is enforced) and refuses one that could cross the cap (`spend.js`
`meteredModel`, a decorator over the model port: it reserves when the request is prepared, so a
refused request is never sent); a JSON-lines ledger keeps the total across restarts. The client sees no usage report, so the charge is an upper bound: the
Copilot's conservative token estimate (2.5 characters per token of text), with every prompt token
at the dearest prompt rate. Prices (USD per million tokens, 2026-10-14; the cache read fixed on
2026-10-24): input $2, output $10, cache hits $0.20, cache writes $2.50, so prompts are charged
at $2.50. The charge counts what streamed back, so thinking billed as output but returned with
no text (Sonnet 5.5 thinks adaptively by default) is not in it: each request's reservation is
bounded by `max_tokens`, the ledger's charge for it is not. `RMC_EVAL_SONNET_MODEL`,
`RMC_EVAL_PRICE_IN` and `RMC_EVAL_PRICE_OUT` override the defaults; `ANTHROPIC_API_KEY` is required.

## Arm A results (no model; 2026-10-14, commit b4fe09f)

`npx vite-node validation/run.js -- arm-a --split dev` (and `--split test --final`). Records:
`validation/baseline/arm-a-dev.json` and, sealed since the holdout gate,
`validation/holdout/arm-a-test.json` (re-recorded after batch 1 at commit 5780c50, below;
the 2026-10-14 records are in git at b4fe09f; the test records again at 4ab9eab for protocol
2026-10-19.2, with the same scores; the dev A+ records again at 3f90a80 for the engine verdicts,
their cluster and modulation readings only, with the same scores and claims; and at 63569f3 for
candidate 8, A+ reading the app's rules, with the same scores and claims apart from
dev-offcentre's δ, 0.201 → 0.202 Å). The readings in the last column are the dev cases'; the test
cases' are in the sealed records.

| Kind | A dev / test | A+ dev / test | Where the default analyses mislead |
| --- | --- | --- | --- |
| offcentre | 1 / 1 | 1 / 1 | They don't: split site found, direction within 7°, δ 0.201 Å (dev) |
| zone_boundary | 0.67 / 0.67 | 1 / 1 | **Ordered antiphase order read as random off-centring** (the PCA cloud is split exactly as for disorder); q missed |
| breathing | 0.50 / 0.50 | 1 / 0.50 | **"No order"**: no page sees cluster breathing. A+ test: breathing found, plus a spurious "shape distortion" from the gyration-axis test (L2) |
| jahn_teller | 0.92 / 0.92 | 0.92 / 0.92 | Split shell seen (Δ 0.31 Å for 0.32 on dev), but not the long-axis family |
| cell_wrap | 0.83 / 0.83 | 1 / 1 | No wrap artefact (the engines handle x = 0), but the rigid dimers are invisible without the cluster engine |
| thermal_null | 1 / 1 | 1 / 1 | They don't over-claim (with the calibrated superstructure test; before it, A+ found "order" in 6 of 12 configurations without any) |
| symmetry_lowering | 0.33 / 0.33 | 0.83 / 0.83 | **The 0.2 Å card says Pm-3m:** wrong group, called centrosymmetric and non-polar. A+ reads the noise floor and gets the group, but has no rule naming the shift |
| sites_only (achiral) | 1 / 1 | 1 / 1 | They don't: the point-group facts are deterministic |
| sites_only (chiral) | 1 / 1 | 1 / 1 | They don't |
| **mean** | **0.81 / 0.81** | **0.97 / 0.92** | |

**The bar.**
- **Arm B must beat A** where the pages mislead: ordered vs disordered splits, collective cluster
  motion, symmetry below the card's tolerance, the JT axis.
- **B must beat A+ to show the language model itself helps.** The engines alone already score
  0.92 on the test set.
- **B must match A on the symmetry-fact cases**, where the card is already perfect and a model
  can only lose. The chirality trap measures exactly that.

The A/A+ rules were set on the dev cases and frozen before the test set was scored. The test-set
A+ spurious distortion was seen and left unchanged.

**Batch 1 (2026-10-17): A+ follows the engines.** A+ is every Copilot engine, so B − A+ stays what
the model adds only if A+ gets the batch's new engines too. The sweep now reads the modulation
tool's facts (zone-boundary order, and order at a general wavevector, review L4) and, when the group
at the noise floor is not the card's, `measure_symmetry_lowering` from the card's 0.2 Å group: each
element shifted above the noise floor is a `symmetry_lowering` claim with its direction and size.
The cluster readings use the engine's shape null (review L2) through the tool. On the dev cases
only dev-symmetry-lowering moves: A+ 0.83 → 1 (Ta 0.065 Å along [0 0 1] for 0.07 Å). The rules
were set on the dev cases; the test set was re-run once batch 1 was committed (`arm-a --split test
--final`, commit 5780c50, protocol 2026-10-17.1), reported and not tuned:

| Kind | A dev / test | A+ dev / test (batch 1) | What moved |
| --- | --- | --- | --- |
| offcentre | 1 / 1 | 1 / 1 | — |
| zone_boundary | 0.67 / 0.67 | 1 / 1 | — |
| breathing | 0.50 / 0.50 | 1 / **1** | test: the spurious "shape distortion" is gone (L2's shape null) |
| jahn_teller | 0.92 / 0.92 | 0.92 / 0.92 | — (the long-axis family is still not named) |
| cell_wrap | 0.83 / 0.83 | 1 / 1 | — |
| thermal_null | 1 / 1 | 1 / 1 | — |
| symmetry_lowering | 0.33 / 0.33 | **1** / **1** | the shift named (L4): Ta 0.065 Å [0 0 1] for 0.07 on dev; the test case's likewise, within the tolerances |
| sites_only (achiral) | 1 / 1 | 1 / 1 | — |
| sites_only (chiral) | 1 / 1 | 1 / 1 | — |
| **mean** | **0.81 / 0.81** | **0.99 / 0.99** | |

**The bar after batch 1.** The engines alone now score 0.99 on both sets, so B cannot beat A+ by
more than 0.01 on the case score: on these cases the question "does the model add to the engines"
reduces to whether B *matches* A+ — reaches the engines' conclusions, no over-claim — and to what
the score does not measure (the answer checks, the evidence gate, number provenance,
reproducibility). B − A (0.81) is still where the Copilot's value over the default pages shows.
The report shows where B matches, adds or loses beside the case score: the feature-only score,
B against A+ component by component and every disagreement with the truth
([The report beside the case score](#the-report-beside-the-case-score-plan-item-16-part-i); on
dev the engines' feature-only score is 0.92, not 0.99).

The dev-only fit-data case (protocol 2026-10-19.6) is not in these tables, whose means are over the
nine paired kinds: A and A+ both score 1 on it, so over the ten dev cases outside the judgment
subset the means are A 0.825 and A+ 0.99. The judgment subset (protocol 2026-10-19.9) is not in
these tables either. Its results are in its PR, the dev records and the `arms.test.js` pins, because
the author of its held-out matches may read this document.

## What the model sees (snapshot tests)

A run means something only while the model is sent the same thing. The tests in
`web_app/frontend/src/__tests__/modelVisible/` write it to readable files in its `__snapshots__/`
folder, with no model and no network:

- `system-prompt/<variant>.txt`: every system prompt `copilotSystemPrompt` builds (tools off; tools
  on with or without the view tools, `ask_user` and the skills' index; and with method `'off'`,
  arm C's, `*.method-off.txt`: tools off, and tools on with or without the view tools and `ask_user`), with a line break
  after each sentence for reading, and the exact length and hash of the string sent;
- `dock/first-request.txt` and `dock/tools.json`: the dock's first request on the Demo run
  (CopilotChat → useCopilot, captured at `streamCopilotChat`): its fields, its messages and the
  `tools` field exactly as sent;
- `harness/first-request.arm-B.txt`, `harness/first-request.arm-C.txt`,
  `harness/tools.arm-B.json` and `harness/tools.arm-C.json`: the same for `runLiveCase` on
  `dev-offcentre`, captured at its model port (the files' header still calls it `callModel`, its
  name before the turn module); a test also checks arm C's request for every skill name, the
  skills' step headings, the prompt's method sentences, `EVIDENCE_RULES` and every piece of tool
  advice; `harness/gate-send-back.arm-B.txt` and `.arm-C.txt`: the evidence gate's send-back on
  the same case (a first answer that concludes "no order" with no analysis run), the message each
  arm's model receives, beside the rule and the missing analyses the record keeps;
- `tool-results/<tool>.json`: one result of every run tool and engine tool on the Demo, plus an
  `ask_user` answer, a `read_skill` and a window that cuts the Ta–Se shell, as the `tool` message
  the agent loop appends (budget applied, checks, `ref`, note); `*.method-off.json`: the Demo's
  results whose checks carry method advice (`find_gr_peaks` on Ta–Se and Ta–Ta, `compute_bond_angles`
  on the whole shell and in the refused cut window, `detect_symmetry`, `analyze_modulations` with
  and without an element named), as arm C reads them;
- `tool-messages/<path>.txt`: every tool message no other file covers, as the raw string the model
  reads, byte for byte and in order: a call that did not run (an unknown tool, a fifth call), the
  questions (refused, malformed, unasked, skipped with and without a default, typed), arguments
  that do not parse or that the schema refuses, a throw, the time limit, a repeat, the Workbench's
  notes (asked, settled, settled then repeated, before the turn), a result fitted to its budget or
  cut, earlier results compacted, and `highlight_cards`' three answers.

The validation's scripted runs (`harness.test.js`, `fitCase.test.js`, `judgment.test.js`,
`turnEndings.test.jsx`) pin their records as files too, in
`web_app/frontend/validation/__tests__/__snapshots__/records/` (`recordFiles.js`: the clock, the
current protocol and each request as sent, by length and hash, as placeholders; the keys the tool
transcript added left out), so a change to what a record holds shows as a diff.

`turnRequests.test.js` (no snapshot file of its own) runs the turn module directly with a
dock-shaped spec (the Demo run, the dock's tools, `providerModel`) and a harness-shaped one
(`dev-offcentre` through the run module, arms B and C): its first requests are the requests in
these files, so the dock and the harness send what the one module builds.

Beside them, `caseBlind.test.js` (no snapshot file) checks that the harness tells the model nothing
about the case:
- the first request on every case, dev and test (the test cases where the loader opens them: in
  CI), in arms B and C, names no case id in any spelling (dashes, underscores, spaces, any case);
- no case kind, and no multi-word case id without its split ("cell-wrap"), appears in that
  request's run-specific parts: the context, the question past the fixed case question, the request
  fields, and every tool string that is not the same in every case. The fixed text (the system
  prompt, the verdict format, the tool names and descriptions) names kinds of order on purpose;
- with every tool the harness offers called once on each dev case (every section, pair, element
  and skill it lists), no request of the turn and no tool result names the case id.

It reads only each case's id and kind, never a test case's parameters or truth (it builds the
sealed cases through the loader, and only in CI).

A change to any of these files is a change to what the model sees: say so in the PR, and bump
`PROTOCOL_VERSION` when it changes what a run means. Refresh them with
`npx vitest run src/__tests__/modelVisible -u` and read the diff before committing; CI never writes
them, so a stale file fails there.

The paths a run's answer quotes by reference (`{{r3.sites.0.shift_A}}`) are pinned too: every live
tool declares the shape of its data on its spec (`dataShape`, never sent to the model), and
`src/__tests__/toolDataShapes.test.js` checks every tool's results on the Demo and on each dev case
against it — no undeclared key, no required key missing, every declared key returned by some call.
A change to a tool's output changes its `dataShape` in the same PR, so it shows in the diff
([llm README](../../web_app/frontend/src/llm/README.md), *The tool contract*).

## Running the live arms (later, on the Mac)

From `web_app/frontend`, with the model loaded in LM Studio (with a context of 48k or more:
`lms load <model> --context-length 49152`) or Ollama. Each run starts with the pre-run check
([The run policy](#the-run-policy-plan-item-12)):

```bash
# smoke pass: one run per dev case
npx vite-node validation/run.js -- live --arm B --model qwen3.8-27b --split dev --runs 1
# dev, both arms, 3 runs (repeat for gemma4:26b and meta/muse-glimmer)
npx vite-node validation/run.js -- live --arm B --model qwen3.8-27b --split dev --runs 3
npx vite-node validation/run.js -- live --arm C --model qwen3.8-27b --split dev --runs 3
# held-out test: from a clean committed tree
npx vite-node validation/run.js -- live --arm B --model qwen3.8-27b --split test --runs 3 --final
# Claude Sonnet 5.5: test set, arm B, $20 cap
ANTHROPIC_API_KEY=… npx vite-node validation/run.js -- live --arm B --model claude-sonnet-5.5 --split test --runs 3 --final
# dev check only: the same with the PCA table in the context (own folder, own report column)
npx vite-node validation/run.js -- live --arm B --model qwen3.8-27b --split dev --runs 1 --pca-table
# side by side
npx vite-node validation/run.js -- report --in validation/runs/<dir> --in validation/runs/<dir> --out validation/runs/report.md
```

A rerun keeps runs already recorded in its folder, so an interrupted batch resumes. Merge the
local-model hardening (`fix/local-model-robustness`) first: its deadlines and loop guards keep a
stalled run from holding the batch.

**Go/no-go before the live runs.**
- **Each LM Studio model at 48k.** Whether the Mac can hold each model with a 48k context (its
  memory, and the model's own maximum context) is checked at the go/no-go, one model at a time:
  load it with `lms load <model> --context-length 49152`, run the smoke pass's first case, and read
  the window in its record. The pre-run check refuses anything less, so a model the Mac cannot hold
  at 48k does not run under this protocol.
- **The final round on LM Studio.** In the smoke records, a run whose last request came back as
  tool calls only ends on the loop's empty-final answer. On LM Studio that would mean
  `tool_choice: 'none'` is not honoured. On Ollama the last request carries no tools (protocol
  2026-10-19.7); check in its smoke records that the final reply reads the earlier calls and
  results from the messages alone and still answers (not verified with a live model).

**Time.** Each run is bounded by the agent loop: 6 rounds, 10 minutes plus the final reply. Since
protocol 2026-10-17.2 the evidence gate can add up to two answer requests to a run (an answer sent
back, then the analysis it asked for inside the same round budget). On
qwen3.8-27b (1–5 minutes a question in the 2026-10-10 spot check; a reasoning round of a few
thousand tokens on a 8–15k-token prompt), expect about 8–12 minutes a run, a few more when the gate
sends an answer back:
- one split, one arm, 3 runs: 27 runs, about 4–5 hours;
- dev and test, arms B and C: about 16–20 hours per local model;
- smoke pass: about 1.5 hours.

**Tokens and cost.** For Sonnet, about 70–90k prompt tokens and 3–6k reply tokens per run, so the
27 test runs of arm B come to about 2.4M prompt + 0.15M reply tokens. That is about $6.30 at list
prices ($4.80 prompt + $1.50 reply), less with cache hits. The tracker's upper bound is $7.50, inside the $20 cap.
