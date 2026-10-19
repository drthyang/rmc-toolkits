# AI Copilot validation: does the Copilot help?

The Copilot's own tests replay scripted answers. This protocol asks a different question: given a
configuration whose local structure is **known**, does the Copilot reach the right conclusion more
often than the Workbench without it, and does it avoid claiming order that is not there? It was
fixed in code (`web_app/frontend/validation/`) before any model run. The tool and skill review
that went with it is [copilot-review.md](copilot-review.md).

## Arms

| Arm | What answers | Code |
| --- | --- | --- |
| **A** | The default Workbench, no agent: the Detected SG card (0.2 Å), the PCA Ellipsoid table, the Displacement Directions page (ν = 10, smoothing 2) and the partial g(r), read by fixed rules (`ARM_A_RULES`). The cluster and superstructure engines have no page, so A cannot use them | `armA.js` `runArmA` |
| **A+** (control) | Every Copilot engine run exhaustively with its own thresholds, joined by fixed rules (`SWEEP_RULES`); no model. **B − A+ is what the language model adds**; B − A mixes that with the new engines | `armA.js` `runArmSweep` |
| **B** | The Copilot as the app assembles it: tools, skills, guardrails (refusals, the Workbench's window questions, answer checks) | `harness.js` |
| **C** | B without the skills (no `read_skill`, no skills index): the ablation | `harness.js` |
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
questions, the evidence gate, the answer checks and the references from the run props. The
harness (`harness.js`) is its second adapter and passes only what is its own: the case's run, its
compute-only tools (the tools over the run props with no cards, and the run's engine tools), the
simulated user, the arm's skills (arm B) and a model port. `run.js` builds that port as the app
does, with `providerModel` (`src/llm/provider/model.js`: the context guard, the warmer retry), and
meters it for an API model (`spend.js` `meteredModel`).

## Cases

Nine kinds, each a **dev** case (for improving tools and skills) and a held-out **test** case (for
reported numbers only): the same physics with other elements, directions, magnitudes or
structures. Each case is an `.rmc6f` supercell (reference sites and cell indices; coordinates
wrapped into [0, 1) as RMCProfile writes them) with its partial g(r). `run.js write` saves each as
`<out>/<case>/rmc_run/`, a run folder the Workbench opens, named as the harness names every run
(the Workbench names a run after its folder), with the truth beside it in `<out>/<case>/truth.json`.

| Kind | Planted (dev / test) | The trap |
| --- | --- | --- |
| offcentre | Two-site disorder: Nb ±0.20 Å along [001] in KNbO₃ / Ge ±0.16 Å along [110] | A split site looks like large thermal motion in the ellipsoid alone |
| zone_boundary | Ordered antiphase Ti ±0.10 Å along c, q = (0,0,½) / Zr ±0.09 Å along a, q = (½,½,½); plus a stronger box-scale wave of the A cation | Identical single-site statistics to two-site disorder; the strongest mode is the decoy |
| breathing | Mo₃ triangles / Fe₄ tetrahedra (Cl anion, edge 2.75 Å; since protocol 2026-10-19.2) whose size fluctuates collectively (sd 4 % / 3.5 %), shape regular | Size is not shape |
| jahn_teller | Isolated CuF₆ (4 × 1.93 + 2 × 2.25 Å) / MnCl₆ (2.48 / 2.85 Å), long axis random over x, y, z | The first peak holds 4 neighbours; the average octahedron is regular |
| cell_wrap | Rigid V₂ / Ti₂ dimers whose first site sits at fractional 0 (half its atoms written near 1) | The cell wrap |
| thermal_null | BaTiO₃ / PbTiO₃ composition, thermal motion only (O flattened across its bond) | Literature off-centring; anisotropic but Gaussian O |
| symmetry_lowering | Every Ta shifted 0.07 Å along [001] (P4mm) / Ti 0.06 Å along [111] (R3m) | The 0.2 Å card still says Pm-3m |
| sites_only (achiral) | Zincblende CuBr, F-43m / sulvanite Cu₃VS₄, P-43m | No inversion does not make a structure chiral |
| sites_only (chiral) | B20 FeSi / MnSi, P2₁3 | A genuinely chiral (Sohncke) group |

Thermal rms: 0.05–0.09 Å. Supercells: 6–8 cells a side, 1,024–2,560 atoms.

**The split is enforced in code** (`__tests__/split.test.js`):
- the test set is pinned by a sha256 of its configurations, partials and truths;
- dev and test share no id, and no tool, skill, app source or doc names a test case, its
  material, its planted lengths or its directions ([the holdout scan](#the-holdout-scan));
- `--split test` needs `--final` and a clean git tree, and every record names its commit;
- a supervisor finding on a test case may only report a bug, with the user's approval.

## The holdout scan

`__tests__/split.test.js` (with `__tests__/holdoutScan.js`) keeps the test cases out of what the
tools, the skills and the docs say. It reads each test case's parameters from `cases.js` at run
time, so no test-case value is written in the test, and a failure prints the line with what matched
masked (■).

- **Tokens, from the parameters:**
  - the id, anywhere;
  - the composition (two or more elements): a formula with exactly that element set, with or
    without subscripts, or, for three or more elements, an element-dash chain naming each once;
  - the cluster written as a formula of its own (a tetrahedron's M₄, a dimer's M₂);
  - a structure-type word that only a test case uses;
  - each planted length in Å (displacement, bond lengths with their split and mean, cluster
    edge, the dimers' shared shift, the lattice parameter) within six words of its element;
  - each planted direction family ([uvw] or <uvw>, any order and sign) and wavevector
    (q ≡ −q ≡ q + G, any order) within six words of its element.

  Lengths, directions and wavevectors count only where they differ from the dev cases of the same
  kind, and thermal widths never count: those are the kind's physics, not the test case.
- **Not hits:** a bare element symbol (they are everywhere); a formula with more or fewer
  elements; a chain that repeats an element (A–X–A is a bond angle, which every study of that pair
  writes); a number with more digits (0.371 is not 0.37) or more than six words from its element.
- **Where:** `web_app/frontend/src` (js, jsx, md, and the snapshots of what the model sees), every
  file in `.claude/skills`, `rmc_toolkits` (py), `docs/` (md) and the top-level md files. Not
  scanned: the validation itself (`web_app/frontend/validation/`) and this document, whose Cases
  table names every test case. [AGENT_PLAN.md](../AGENT_PLAN.md) is scanned; it pointed at the
  open leak below without naming it.
- **Checks on the scanner:** every test case yields more than its id, and the scan finds every
  test case in its own truth (summary and traps, prose written from the same parameters). Its own
  rules are pinned on made-up cases (`holdoutScan.test.js`).

Every hit must sit in one of two lists in `split.test.js`, and each entry pins its exact number of
hits, so a new mention, or a removed one, fails until the entry is reviewed. Any hit in neither list
fails the test.

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
module (`src` outside `__tests__`) names it, so none can import it; and that no snapshot contains
its matched text.

| File | Case · token | Hits | Why it is not model-visible |
| --- | --- | --- | --- |
| `web_app/frontend/src/__tests__/fixtures/symmetryStructures.js` (65) | symmetry_lowering · composition | 1 | Test only: a symmetry-finder fixture of a real orthorhombic perovskite with the case's composition, imported only by the symmetry tests; unrelated to the planted shift |

It was found on 2026-10-19 and reviewed against that rule. Four entries found that day for the test
breathing case's old cluster (the llm README, a triplets engine test, ai-assistant.md and
bond-geometry.md, all doc or test only) left the list with protocol 2026-10-19.2: they name the
skills' reference compound, which is no longer a test case.

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
- unsafe actions (a window marked as the user's choice when no question was asked; must be 0);
- unknown-tool and failed calls;
- **number provenance:** each number with a decimal point or a unit, matched to a value the
  model was shown, or "unsourced";
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

**Protocol 2026-10-17.2 (THY-18).** The live arms run with the evidence gate, numbers by reference
and the computed verdict, exactly as the app does. The case score is unchanged; the three are
reported beside it (`report.js`: one more table per split) and in each transcript (the answers the
gate sent back, the answer as written and as filled, the computed outcome with its reasons).

**Protocol 2026-10-19.1 (the case-id leak).** Until this version the harness named each run, and
its `.rmc6f` and partials files, after the case id. The run name reaches the model as the
context's `run`, in the first message, so every arm B and C run was told its case kind
(`"run": "dev-zone-boundary"`). Runs are now named `rmc_run`, the same for every case. Live runs
under an earlier protocol, if any, are not comparable with later ones and must be run again. Arms
A, A+ and D have no model and are unaffected; their records keep protocol 2026-10-17.1.

**Protocol 2026-10-19.2 (test-breathing recomposed).** The held-out breathing case planted Nb₄
tetrahedra with a Se anion: the skills' reference compound, whose published shells
`.claude/skills/coordination-shells/lacunar-spinel-shells.md` gives and `read_skill` serves to arm
B. That leaked the case's kind to arm B only. The fix goes in the test case; the skill is the
user's text and is unchanged.
- **The new composition:** Fe₄ tetrahedra with a Cl anion, a = 6.90 Å, edge 2.75 Å (Fe–Fe
  ≈ 2.73–2.78 Å in Fe₄S₄ cubane cores). Shape, breathing (sd 0.035), thermal sd (0.05 Å) and
  supercell (8³) are unchanged, and a/edge = 2.509 keeps the old shell separation. The shortest
  Fe–Fe between clusters is 5.32 Å (1.94 × edge). A fresh author chose it, without seeing dev
  results, dev parameters or arm records. No Fe, Cl, Fe–Cl, Fe₄, iron or chlorine word appears in
  any skill (main or PR #71's drafts) or in the model-visible snapshots.
- **The seed, by a rule fixed before looking:** seeds 731, 732, … in order; the first whose drawn
  ε sample sd over the 512 clusters (the generator's own `cellState` draws) lies in
  [0.033, 0.037]; at most 20 tries. 731 gave 0.03273 (out) and 732 gave 0.03696 (in), so the seed
  is 732, used by no other case. The rule depends only on the case's own draws, not on any model or
  arm result.
- **Difficulty, honestly:** the shorter edge makes the planted radius spread 7 % smaller
  (0.035 R = 0.059 vs 0.063 Å). The realised spread goes the other way: seed 732's ε lie above the
  planted sd (0.03696) and the old seed's below (0.03362), so the realised radius spread is
  0.0622 Å against the old case's 0.0607 Å (2.5 % larger). A+'s cluster engine reads the breathing
  at z = 30, against 22 before.
- **Approved by the manager on the user's delegation, 2026-10-19**
  (`validation/findings/decisions.jsonl`). The user may reverse it: keep their compound in the test
  case and disclose it as a known arm-B advantage.
- **Nothing is lost:** no arm B or C run has touched the test set. `TEST_SET_SHA256` changes, so
  any result on the old test set would not be comparable.
- **Arms A and A+ on the test set**, re-recorded with `--final` at 4ab9eab: only the three
  test-breathing records change (readings, the verdict's element, the truth), and no score moves.
  A scores 0.50 (no page sees cluster breathing). A+ scores 1: Fe breathing found (breathing sd
  0.107 vs 0.0635 Å for independent atoms, z = 30), no shape distortion (edge spread 0.152 vs
  0.154 Å). D scores 1. The means stay A 0.81 and A+ 0.99. The dev records are unchanged: a re-run
  equals `arm-a-dev.json` record for record, apart from protocol and commit.
- **The dev set overlaps the skills, and may.** dev-breathing (Mo₃ triangles, O anion) was checked
  against the same rule. No Mo₃, Mo–O or Mo appears in the skills on main or in the snapshots.
  PR #71's draft `local-order/order-signatures.md:25` reads "Breathing Mo₃ triangles (size sd
  4 %)": dev-breathing's cluster and its planted breathing. The dev set is for development, and
  the skills learn from it. Once that draft merges, read an arm B result on dev-breathing with this
  in mind; the reported numbers come from the test set.
- **The holdout scan's open-leak list is now empty**, and the four allow-list entries for the old
  cluster have left it ([the holdout scan](#the-holdout-scan)).

**Protocol 2026-10-19.3 (the turn module).** The harness had assembled each question itself,
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
   the PCA table), so no default record moves.
4. **Every request as sent.** The record keeps `requests` (below).

Changes 1 and 2 change what a run sends, so the protocol moves. Live runs under an earlier
protocol that looped or met a model without tools are not comparable; arms A, A+ and D have no
model and are unaffected.

## Records, transcripts, findings

Each live run writes:
- `<case>__<arm>__<model>__r<n>.json`: the view the run was opened with and what the model was
  told about it (below), every step with who decided it (the model or the Workbench), the
  reasoning and text the model wrote before the call, the arguments, the checks and the result it
  read, then the final answer, verdict, score, answer checks, provenance, safety, and the truth;
- `.md`: the same as a self-contained transcript, with numbered steps and the planted truth last.

**What the run had (plan item 1.3).** The app's Copilot can be handed more than a case gives the
harness, so each record states it:
- `view`: the run module's normalised view, `{ toleranceA, hiddenPlots, pcaSites, liveData,
  ladder }`. By default it is `HARNESS_VIEW`: 0.2 Å with the ladder, nothing hidden, no PCA table,
  Live Data off. In the app the tolerance is the user's ladder pick, and the PCA table joins once
  the PCA Ellipsoid page has computed it;
- `context`: `sections`, the context sections the model read (`harness.js` `CONTEXT_SECTIONS`, in
  `buildRunContext`'s order); `absent`, the sections it lacked; `earlier_turns: 0`, since every run
  is one fresh question (the app sends up to eight earlier turns). A case has `structure`,
  `symmetry`, `pair_correlations` and `datasets` (its partial g(r)). It has no
  `configuration_optimization` (no move counters), `run_settings`, `partial_weights` or
  `convergence`, and no `pca_displacements` unless the PCA table is on. Read from the context as
  the model was shown it (the turn's `start`);
- `requests` (since protocol 2026-10-19.3): every model request as it was sent, in order, as
  `{ round, why, sent }`: `why` is `first`, `tools` (the results of the round before), `reminder`
  (after an empty or looping reply), `gate` (an answer the evidence gate sent back), `final` (the
  analysis budget used: `tool_choice` 'none') or `no_tools` (the retry of a model that refuses
  tools); `sent` is the request as the port prepared it (`providerModel`: base URL, model id,
  temperature, messages, tools, tool choice; never the API key). A request the spend cap refused
  was never sent and is not there. One the context guard refuses (it would not fit the model's
  window) is there, followed by the run's `error`: the guard decides as the request goes out, not
  when it is prepared (the model-port PR, candidate 3, moves that decision into `prepare`).

The transcript prints both under its header. **The PCA-table option** (`run.js live --pca-table`,
dev set only; `runLiveCase({ pcaTable: true })`) gives each run the PCA table from its own engine,
the rows `list_sites` reads, as after a visit to the PCA Ellipsoid page. It is a dev check of what
that section changes and is off by default. Such runs write to their own folder (`…_pca-table`),
and `report` gives them their own column (`<model> · <arm> + PCA table`). It changes the context the
model reads, so those runs are not comparable with the default runs.

`report` puts every record side by side (case × model·arm), keeps failed runs as score 0, and
lists every run used.

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
L1–L5 (L4 is two analyses), six findings; the five skill proposals are added by the skills PR,
which waits for the user's review. `split.test.js` validates the file.
Decisions on the validation itself take the same form in `validation/findings/decisions.jsonl`,
also validated by `split.test.js`: so far one, the recomposed test breathing case (protocol
2026-10-19.2), approved by the manager on the user's delegation on 2026-10-19. The schema's
approver is the user, so the delegation is stated in the entry's evidence and proposed fix.

## Models and spend

Each model gets the same cases, prompt, tools and scoring, and the results are reported side by side:
- qwen3.8-27b (LM Studio, 32k);
- gemma4:26b (Ollama);
- meta/muse-glimmer (LM Studio, 28B, tool use, up to 128k);
- Claude Sonnet 5.5 (API).

The API model (`claude-sonnet-5-5`) runs the test set only, arm B, 3 runs per case, under RMC's
$12 share of the cycle's $40 cap. The spend tracker prices each request before sending it (prompt
estimate + the largest reply) and refuses one that could cross the cap (`spend.js`
`meteredModel`, a decorator over the model port: it reserves when the request is prepared, so a
refused request is never sent); a JSON-lines ledger keeps the total across restarts. The client sees no usage report, so the charge is an upper bound: the
Copilot's conservative token estimate (2.5 characters per token of text), with every prompt token
at the dearest prompt rate. Prices (USD per million tokens, 2026-10-14): input $2, output $10,
cache hits $0.10, cache writes $2.50, so prompts are charged at $2.50. `RMC_EVAL_SONNET_MODEL`,
`RMC_EVAL_PRICE_IN` and `RMC_EVAL_PRICE_OUT` override the defaults; `ANTHROPIC_API_KEY` is required.

## Arm A results (no model; 2026-10-14, commit b4fe09f)

`npx vite-node validation/run.js -- arm-a --split dev` (and `--split test --final`). Records:
`validation/baseline/arm-a-{dev,test}.json` (re-recorded after batch 1 at commit 5780c50, below;
the 2026-10-14 records are in git at b4fe09f; the test records again at 4ab9eab for protocol
2026-10-19.2, with the same scores).

| Kind | A dev / test | A+ dev / test | Where the default analyses mislead |
| --- | --- | --- | --- |
| offcentre | 1 / 1 | 1 / 1 | They don't: split site found, direction within 7° / 10°, δ 0.201 / 0.162 Å |
| zone_boundary | 0.67 / 0.67 | 1 / 1 | **Ordered antiphase order read as random off-centring** (the PCA cloud is split exactly as for disorder); q missed |
| breathing | 0.50 / 0.50 | 1 / 0.50 | **"No order"**: no page sees cluster breathing. A+ test: breathing found, plus a spurious "shape distortion" from the gyration-axis test (L2) |
| jahn_teller | 0.92 / 0.92 | 0.92 / 0.92 | Split shell seen (Δ 0.31 / 0.38 Å for 0.32 / 0.37), but not the long-axis family |
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
| symmetry_lowering | 0.33 / 0.33 | **1** / **1** | the shift named (L4): Ta 0.065 Å [0 0 1] for 0.07; Ti 0.055 Å [1 1 1] for 0.06 |
| sites_only (achiral) | 1 / 1 | 1 / 1 | — |
| sites_only (chiral) | 1 / 1 | 1 / 1 | — |
| **mean** | **0.81 / 0.81** | **0.99 / 0.99** | |

**The bar after batch 1.** The engines alone now score 0.99 on both sets, so B cannot beat A+ by
more than 0.01 on the case score: on these cases the question "does the model add to the engines"
reduces to whether B *matches* A+ — reaches the engines' conclusions, no over-claim — and to what
the score does not measure (the answer checks, the evidence gate, number provenance,
reproducibility). B − A (0.81) is still where the Copilot's value over the default pages shows.

## What the model sees (snapshot tests)

A run means something only while the model is sent the same thing. The tests in
`web_app/frontend/src/__tests__/modelVisible/` write it to readable files in its `__snapshots__/`
folder, with no model and no network:

- `system-prompt/<variant>.txt`: every system prompt `copilotSystemPrompt` builds (tools off; tools
  on with or without the view tools, `ask_user` and the skills' index), with a line break after each
  sentence for reading, and the exact length and hash of the string sent;
- `dock/first-request.txt` and `dock/tools.json`: the dock's first request on the Demo run
  (CopilotChat → useCopilot, captured at `streamCopilotChat`): its fields, its messages and the
  `tools` field exactly as sent;
- `harness/first-request.arm-B.txt`, `harness/first-request.arm-C.txt` and
  `harness/tools.arm-B.json`: the same for `runLiveCase` on `dev-offcentre`, captured at its model
  port (the files' header still calls it `callModel`, its name before the turn module);
- `tool-results/<tool>.json`: one result of every run tool and engine tool on the Demo, plus an
  `ask_user` answer, a `read_skill` and a window that cuts the Ta–Se shell, as the `tool` message
  the agent loop appends (budget applied, checks, `ref`, note).

`turnRequests.test.js` (no snapshot file of its own) runs the turn module directly with a
dock-shaped spec (the Demo run, the dock's tools, `providerModel`) and a harness-shaped one
(`dev-offcentre` through the run module, arms B and C): its first requests are the requests in
these files, so the dock and the harness send what the one module builds.

Beside them, `caseBlind.test.js` (no snapshot file) checks that the harness tells the model nothing
about the case:
- the first request on every case, dev and test, in arms B and C, names no case id in any spelling
  (dashes, underscores, spaces, any case);
- no case kind, and no multi-word case id without its split ("cell-wrap"), appears in that
  request's run-specific parts: the context, the question past the fixed case question, the request
  fields, and every tool string that is not the same in every case. The fixed text (the system
  prompt, the verdict format, the tool names and descriptions) names kinds of order on purpose;
- with every tool the harness offers called once on each dev case (every section, pair, element
  and skill it lists), no request of the turn and no tool result names the case id.

It reads only each case's id and kind, never a test case's parameters or truth.

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

From `web_app/frontend`, with the model loaded in LM Studio or Ollama:

```bash
# smoke pass: one run per dev case
npx vite-node validation/run.js -- live --arm B --model qwen3.8-27b --split dev --runs 1
# dev, both arms, 3 runs (repeat for gemma4:26b and meta/muse-glimmer)
npx vite-node validation/run.js -- live --arm B --model qwen3.8-27b --split dev --runs 3
npx vite-node validation/run.js -- live --arm C --model qwen3.8-27b --split dev --runs 3
# held-out test: from a clean committed tree
npx vite-node validation/run.js -- live --arm B --model qwen3.8-27b --split test --runs 3 --final
# Claude Sonnet 5.5: test set, arm B, $12 cap
ANTHROPIC_API_KEY=… npx vite-node validation/run.js -- live --arm B --model claude-sonnet-5.5 --split test --runs 3 --final
# dev check only: the same with the PCA table in the context (own folder, own report column)
npx vite-node validation/run.js -- live --arm B --model qwen3.8-27b --split dev --runs 1 --pca-table
# side by side
npx vite-node validation/run.js -- report --in validation/runs/<dir> --in validation/runs/<dir> --out validation/runs/report.md
```

A rerun keeps runs already recorded in its folder, so an interrupted batch resumes. Merge the
local-model hardening (`fix/local-model-robustness`) first: its deadlines and loop guards keep a
stalled run from holding the batch.

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
prices ($4.80 prompt + $1.50 reply), less with cache hits. The tracker's upper bound is $7.50, inside the $12 cap.
