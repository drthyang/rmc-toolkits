# Copilot plan: architecture, tools and skills before the live runs

The validation ([copilot-validation.md](algorithms/copilot-validation.md)) measures the AI Copilot
as it ships. Before arms B and C run on the test set, this plan fixes what would make the numbers
measure the harness instead of the Copilot, puts each piece of science in one place, and brings the
skills into line with the tools. Then the code is frozen, arms A and A+ are run again, and the
validation runs as written.

It builds on the [tool and skill review](algorithms/copilot-review.md) (2026-10-14), validation
batch 1 and THY-18 (PRs #69 and #70, merged), and a survey of `web_app/frontend/src` (the Copilot,
its tools and the workers), `web_app/frontend/validation`, `rmc_toolkits` and `.claude/skills` on
`main` at `b5b3820` (2026-10-18). Every `file:line` below is on that commit. Sections B and C
summarise the architecture review of 2026-10-18 and the user's decisions of 2026-10-19. The skill
drafts of batch 1 wait in PR #71 for the user's review.

## Ground rules

- **The test set stays held out.** Every change here is designed and checked on the dev set, the
  scripted replays and CI. Arms B and C do not run on a test case until the freeze (phase 6). Arms
  A and A+ are deterministic and are re-run at the freeze, since several items change the engines
  A+ sweeps.
- **Approval.** A new tool, a skill change, a tool description or prompt rule that restates a
  skill, and a threshold change that alters the science all need the user's approval. Items waiting
  for it are marked **(approval)**.
- **No API spend.** Local models and scripted CI only, until the user says otherwise. The
  Sonnet arm's share ($20 since the user's decision of 2026-10-24, $12 before) stays unspent.
- **Each PR** passes CI (Python on 3.9 floors, 3.11, 3.13; frontend lint, vitest, build), updates
  its docs, and leaves the app working. A refactor PR keeps the model-visible text unchanged: a
  snapshot of the prompt, the tool list and the descriptions (phase 1.7) shows any change in the
  diff.

## Status (plan-status 5, 2026-10-24)

**Scope.** On 2026-10-19 the user widened the work to the whole plan, through phase 6, and
delegated approvals to the manager. Still the user's: skill text (#71, 5.2, 5.4, and 5.5's
on-demand references), live model runs and API calls, and releases. Items approved under that
delegation are marked *approved 2026-10-19* below. Each one keeps the guardrails: a known-answer
test, arms A and A+ re-run on dev with every moved case reported, the snapshot diff shown in its
PR, and a decision record in `validation/findings/decisions/` and here.

**Merged so far:**

| Item | PR | What |
| --- | --- | --- |
| 1.7 | #73 | snapshots of what the model sees |
| 7 | #74 | generated goldens at the orientation and site seams |
| a | #75 | the case-id leak fixed (protocol 2026-10-19.1) |
| 1.5 | #76 | the holdout scan widened |
| 3.3 | #77 | every live tool declares its `dataShape` |
| 2, 1.3 | #78 | the run module; records keep their view |
| — | #79 | the held-out breathing case recomposed (protocol 2026-10-19.2) |
| 5.1 | #80 | the skill lint |
| 5.5 (test) | #81 | the 32k budget test |
| 4.1 | #82 | `list_sites` shows the split (protocol 2026-10-19.3) |
| 3.5 | #83 | a corpus of real answers for the claim regexes (0 false alarms in 4 readings, 1 miss) |
| 1, 1.1 | #84 | the turn module; the harness's warmer retry (protocol 2026-10-19.4) |
| — | #85 | plan status 1 |
| 1.2 | #86 | one run policy, recorded; LM Studio below 48k refused before a run (protocol 2026-10-19.5) |
| 4.6 | #87 | a fit-data dev case (protocol 2026-10-19.6) |
| 3 | #88 | the provider quirks behind the model port (protocol 2026-10-19.7) |
| 5 | #89 | the engines return their own verdicts (clusters, modulations) |
| 1.4 | #90 | state lifecycle: a new configuration ends the earlier answers; the pick resets per run |
| 1.6 (i) | #91 | the report compares B with A+ per component; feature-only scores |
| 9 | #92 | one request contract, checked against every adapter |
| — | #93 | an App-level smoke test of the dock; two Live Data fixes |
| 6 | #94 | one shell reading per pair, per run; the first peak at the inter-peak minimum (protocol 2026-10-19.8) |
| — | #95 | plan status 2 |
| 2.6 | #96 | one PCA fit subsample in both runtimes above the 20,000-copy cap |
| 1.6 (ii) dev | #97 | five judgment dev cases; the spec for the blind held-out author (protocol 2026-10-19.9) |
| 1.8 | #98 | arm C without the method: a structural switch; the gate names no analysis in C (protocol 2026-10-19.10) |
| 8 (2.3, 2.4) | #99 | arm A+ reads the app's rules: the split site, the cluster elements, the lower group's tolerance |
| 1.6 (ii) test | #100 | five held-out judgment cases written blind, with a disclosure (protocol 2026-10-19.11) |
| 4 | #101 | one verification module answers "is this answer verified?" |
| 4.3, F2, F3 | #102 | each ordered mode's polarisation; no value rounded onto its bar; one Šidák family over every element (protocol 2026-10-19.12) |
| 10 | #103 | one turn policy decides what ends a turn |
| — | #104 | the held-out test set sealed behind one loader; no held-out value in any doc; the scan covers docs, decisions and fragments |
| — | #105 | plan status 3 |
| F5 | #106 | a split M₃/M₄ cluster shell does not join: a known limit, with its evidence |
| 11 | #107 | result statuses on the tool contract, the event stream as the transcript; "asked before" and provenance read only what was shown or told (protocol 2026-10-19.13) |
| 4.5 | #108 | the structure's basis from `/api/structure`: the Detected SG card, the CIF and the symmetry tools for server folders |
| 4.2 | #109 | a joined shell's split by number: long − short ± its sampling error and each part's share (protocol 2026-10-19.15) |
| 2.5 | #110 | the Displacement Directions page at the engine's resolution, its readouts from raw counts; a peak direction reported ± its cell's radius (protocol 2026-10-19.17) |
| — | #111 | Flask mode: the Copilot gets the plot data (the same 11 analyses and context as the browser path) |
| 4.4, F1 | #112 | an independent-cell null for the modulation tests; a test for diffuse order between cells (protocol 2026-10-19.16) |
| 6 (cut rule) | #113 | one window-cut rule, the g(r) CN margin (approved by the user; protocol 2026-10-19.18) |
| 1.6 (ii) dev | #115 | the dev marginal situation re-set for the corrected null (protocol 2026-10-19.19) |
| 1.6 (ii) test | #116 | the held-out marginal pair re-set by a fresh blind author (protocol 2026-10-19.20) |
| — | #117 | plan status 4 |
| 5.2 | #71 | the batch-1 skill drafts: the review's proposals and a new local-order skill |
| 5.2, 5.6 | #114 | the skills' facts synced with the tools, and the five interviews: every skill "Interviewed 2026-10-24" |
| — | #119 | the fit-assessment interview's rules in code: convergence on every χ² term and the total, the 5 % constraint, two answer checks (protocol 2026-10-24.1) |
| — | #121 | the interview answers of four skills as checks: coordination shells, local symmetry, displacement anisotropy, local order (protocol 2026-10-24.2) |
| — | #122 | the χ² clamp test compared to 1e-12 relative (the oldest-deps CI job) |
| — | #118 | `polyhedral_distortion`: a polyhedron's distortion indices, per polyhedron and for the average structure (protocol 2026-10-24.3) |
| — | #120 | the G(r) error ribbons: F(Q)'s σ propagated through RMCProfile's transform (protocol 2026-10-24.4) |

**Frozen (2026-10-24)** at commit 226f6f2 (protocol 2026-10-24.4), after the changelog fold
(#123): arms A, A+ and D re-recorded with `--final` from a clean tree on dev and test, no score
moved (`validation/baseline/freeze.json`); no tag (tags are the user's).

**Then:** the live runs (arms B and C), on the user's go/no-go.

**Waiting for the user:**
- the live-run go/no-go, after the freeze;
- the user's post-v1 list:
  - the neighbouring-ν choice for a direction finding;
  - the "most anharmonic" ranking;
  - A5: whether a tool's own refusal (a bond-angle window that cuts the shell) counts as the
    analysis behind a coordination claim;
  - the eval gaps no check can catch.

**Decisions** (the manager, on the user's delegation, 2026-10-19):
- **1.8:** arm C removes the whole method, through a structural switch in the turn module.
  - **Q4:** in C the evidence gate still runs and checks, so the guardrail is the same in both
    arms. Its send-back names no analysis and no tool ("this conclusion needs its evidence from
    an analysis run in this question"), and no other refusal in C teaches the method.
- **#6, the shell catalogue:**
  - The first-peak window moves to the inter-peak minimum (approved). r ± FWHM overcounts into the
    second peak: on a constructed 3+3 split the minimum gives 3.03 and r ± FWHM 3.4–4.0, and the
    Demo's Ta–Se reads 3.25, not 3.86.
  - The window-cut rule: the user approved the CN margin on 2026-10-19 (#113), one rule for the
    Workbench's question and the tool's refusal; the coordination counts decide only where no shell
    can be compared.
- **Findings F1–F5** (from building the judgment cases): F1, the modulation null for sites that share
  a cell-level variable, is fixed in 4.4 and calibrated at or below α; F2 and F3 are done (#102); F4
  did not reproduce (#99); F5: the shell join is left as it is (accepting CN 2 and 3 helps no case
  and breaks one), a documented known limit (#106).
- **The sealed test set (#104):** the test definitions load only through one loader, on the final
  path and in CI. No doc, decision record, fragment or changelog states a held-out value. Builders
  read only the working tree: no history (log -p, show, blame, diff) of the validation docs,
  decisions, fragments, the changelog or the sealed folder. A builder that has seen a held-out value
  does no further science, threshold or scoring work.
- **#11:** "asked before" counts only a question the user saw or answered, and provenance credits
  only what the model was told. No evidence never counts. The freeze report shows any test-split
  number they move beside the old reading.
- **4.4:** approved, provided it is calibrated on thermal-only seeds (the observed false-positive
  rate at or below α), has a known-answer test on a planted spread, and uses dev cases only. Built
  (#112): at or below α on 1,000 fresh configurations of every dev kind, shared-variable kinds
  included.
- **2.5 (#110):** the page defaults to the engine's resolution. The direction resolution is a
  measurement limit: a peak direction is reported ± its cell's circumscribed radius, and above the
  15° tolerance as not resolved finer. Arm A reads the default page as a user would; arm A+ prefers
  PC1 where the peak cell is wider than 15°.
- **The marginal situation:** F1's calibrated null raised the test's power, so the dev and held-out
  marginal pairs were re-set (#115, #116) from the spec's rewritten knob; the held-out pair by a fresh
  blind author, its blindness checked from its tool log. The user may reverse it before the freeze.
- **The user's pre-freeze decisions (2026-10-19):** the skills' facts fixed now, the interviews after
  v1 (#71, #114); the window-cut rule on the CN margin (#113); the live runs approved for after the
  freeze — local models first, then Sonnet 5.5 within the $40 cap, one heavy model at a time at 48k
  or more, synthetic cases only, with streaming and the divider after a configuration change the
  first things to watch.
- **The skill interviews (the user, 2026-10-24):** all five skills are now the user's own answers
  (#114), and their rules are in code: every-term convergence and the 5 % constraint (#119), the
  four skills' answer checks (#121), `polyhedral_distortion` for an extended network (#118), and
  the G(r) error ribbons for residuals judged against σ(r) (#120). The records:
  `2026-10-24-fit-assessment-interview.json`, `2026-10-24-coordination-shells-interview.json`,
  `2026-10-24-local-symmetry-interview.json`, `2026-10-24-displacement-anisotropy-interview.json`,
  `2026-10-24-local-order-interview.json`, `2026-10-24-polyhedral-distortion.json` and
  `2026-10-24-ft-error-propagation.json`, in `validation/findings/decisions/`.
- **5.4:** arm B's descriptions keep their method advice until the skill text carries it; arm C
  drops it through 1.8.
- **The held-out breathing case:** recomposed, because its cluster was the skills' reference compound,
  which `read_skill` serves to arm B. The user may reverse this (its record is in
  `validation/findings/decisions/`).
- **The context budget (Q3):** the guard, its estimate (1.4 characters per token for tool
  messages), the 8,192-token reserve, the budgets and the prompt stay unchanged; refusing rather
  than cutting is the behaviour wanted.
  - The live runs load each LM Studio model with at least a 48k window. RunPolicy records the
    minimum per server, the harness checks the loaded window before a run and refuses below it,
    and the record stores the window (Ollama's num_ctx included).
  - Whether the Mac can hold each model at 48k is checked at the live-run go/no-go.
- **Process:**
  - PRs add a changelog fragment in `docs/changes/unreleased/` instead of editing CHANGELOG.md
    (folded into the changelog at phase 6).
  - Plan edits land in plan-status PRs like this one.
  - Decision records are one file each.

## What is wrong now

### A. The validation would partly measure the harness

- **The engine path is the app's; three things around it are not.**
  - **The engine.** The harness runs the engines through the same message handler as the
    Copilot's worker, in-thread (`validation/runFolder.js:28`).
  - **The tools.** The same tool builders (`validation/harness.js:82-87`).
  - **The THY-18 checks.** The evidence gate and the reference fill are built as the app builds
    them (`harness.js:130-131`, `src/llm/useCopilot.js:175-180`).
  - **The context.** It differs from what a user's Copilot sees: `pcaSites`, `runSettings` and
    `rValueFile` are null in the harness (`harness.js:77-78`). In the app the PCA table joins the
    context once the PCA Ellipsoid page has loaded (`src/App.jsx:687,739`), and the Detected SG
    tolerance is the user's ladder pick (`src/components/Dashboard.jsx:423`), not the harness's
    0.2 Å (`runFolder.js:15,24-27`). **Stated (1.3, with the run module):** the harness opens each
    case through the run module's in-process adapter with the app's default view, and every record
    keeps that view (`record.view`), its context sections and `earlier_turns: 0`; `--pca-table`
    adds the PCA table for a dev check.
  - **The prompt and history.** The harness offers no view tools, so the model reads the
    compute-only rules instead of `VIEW_RULES` (`harness.js:90-92`,
    `src/llm/prompts/copilot.js:151`), and it never has earlier turns (`harness.js:91`); the
    app sends up to 8 (`copilot.js:125`, `useCopilot.js:137`). These are known and documented
    choices; the record does not say them.
- **The harness tells the model the case kind.** The run is named after the case id
  (`validation/runFolder.js:30`), and the name reaches the model as `context.run` in the first
  message (`src/llm/context/runContext.js:422`, `src/llm/prompts/copilot.js:109,163`). So every
  arm B and C run is told which kind of case it is ("dev-zone-boundary"). `harness.js:13-14`
  promises the opposite, and the holdout scan reads only source text, so it cannot see this.
  Fixed first (B and C, item a). **Fixed (protocol 2026-10-19.1):** every case is loaded as one
  run, `rmc_run`, with files `rmc_run.rmc6f` and `rmc_run_PDFpartials.csv` (`runFolder.js`
  `RUN_NAME`, `RUN_FILES`), and `run.js write` names the run folder the same way.
  `src/__tests__/modelVisible/caseBlind.test.js` checks the first request of every case (dev and
  test, arms B and C) and every tool's results on the dev cases for the case id and kind.
- **The warmer retry is lost in the harness.** After a looping reply the loop asks again with
  `antiRepeat` (`src/llm/copilot/agentLoop.js:374`).
  - **The app** raises the temperature to at least 0.6 for that request (`useCopilot.js:48,148-151`).
  - **The harness** passes the flag into `streamCopilotChat` (`validation/run.js:159-160`), whose
    signature has no such field (`src/llm/provider/contextGuard.js:220-221`), so the retry runs at
    0.2.
  - **The effect:** a model that loops recovers in the app and not in the validation.
- **Arm C still gets much of the method.**
  - **The ablation itself is structural.** No `read_skill` and no skills index
    (`harness.js:81,86,92`); no skill name in any model-visible tool text (the one mention is a
    comment, `src/copilotTools.js:355`).
  - **But the method travels in the prompt and the tools:**
    - **The prompt.** `EVIDENCE_RULES` states which analysis each conclusion needs, the
      local-order method in prose (`copilot.js:32-39`). `TOOL_RULES` names the tools for a
      distortion (`copilot.js:55-58`).
    - **The tool descriptions** restate coordination-shells: "use the whole first coordination
      shell (find_gr_peaks' first_shell.window_A)" (`copilotTools.js:1165`), "take the window from
      find_gr_peaks" (`copilotTools.js:1242-1243`).
    - **The gate's send-backs** spell out what to run (`src/llm/copilot/evidence.js:162`).
  - **So** C measures the skills' text over a prompt that already carries most of their rules.
    That is a legitimate arm, "skills as written", but not "no method", and the doc does not say
    which it is.
  - **Decided (1.8, 2026-10-19):** arm C removes the whole method, as MATERIA's 1.3 does. With
    skills off, the model sees no skill names, stage lists, `EVIDENCE_RULES`, method sentences in
    `TOOL_RULES` or method advice in tool descriptions. It is a structural switch in the turn
    module (B and C, item C), not text edits.
- **Models run under different rules,** so a model comparison is partly a comparison of these:
  - **Context sizing.** Ollama gets a sized `num_ctx` and a `num_predict` (`contextGuard.js:15-19`,
    `src/llm/provider/ollama.js:82`). LM Studio is refused a turn its loaded window cannot hold
    (`contextGuard.js:20-24`). Cloud servers are sent unchanged, with no reply cap
    (`contextGuard.js:25`).
  - **The limits and the guards are the same for all.** The temperature (0.2, `validation/models.js:53`,
    `src/llm/settings.js:18`), the round, time and character budgets (`agentLoop.js:87`, `DEFAULT_LIMITS`),
    the stall and loop guards, and dropping earlier reasoning (assistant messages carry `content`
    only, `agentLoop.js:487,504-505`).
  - **The final round** sends `tool_choice: 'none'` (`agentLoop.js:372`). No test or run shows
    that the Anthropic endpoint (`models.js:33`) accepts it.
  - **No record states the policy it ran under.** Context size, reply cap or server kind are
    missing from the run record.
- **State outlives what it describes.**
  - **Live Data.** A save bumps `configEpoch` but not the Copilot's `runKey` (`App.jsx:85,104`), so
    the conversation keeps up to 8 earlier answers about the previous configuration
    (`useCopilot.js:57,137,222-227`) while the context describes the new one.
  - **The Detected SG pick** (`App.jsx:92`, set at `src/components/ModelSummary.jsx:356`) is never
    reset when another run opens, so a tolerance chosen for one run sets the next run's context.
  - The validation runs one question per fresh run, so neither shows there. Both shape what users
    get.
- **The holdout scan is narrow.**
  - **What it checked.** Only the test-case ids, only in `src` (js/jsx), `.claude/skills` and
    `rmc_toolkits` (`validation/__tests__/split.test.js` before 1.5).
  - **What it missed.** The test cases' materials, formulas, magnitudes and directions, `docs/`,
    and the llm README.
  - **One overlap is a leak, pending the user:** the test breathing case's cluster is the
    cluster of the skills' reference compound, whose published distortion they describe
    (`.claude/skills/coordination-shells/lacunar-spinel-shells.md:9,31-32`). `read_skill` loads
    that file, so it reaches arm B, and a model may expect a distortion there. The skill text is
    the user's, so the decision is too.
  - **Widened (1.5, 2026-10-19).** The scan derives each test case's id, composition, cluster,
    structure type, planted lengths, directions and wavevector from the test definitions at run time, and
    looks across `src`, every skill file, `rmc_toolkits`, `docs/` and the top-level docs. Only
    text no model can read may be allow-listed (five test- or doc-only hits found that day). The
    leak above is listed apart, printed on every run and pinned line by line; both lists are in
    [copilot-validation.md](algorithms/copilot-validation.md#the-holdout-scan). Since #104 the
    definitions are sealed (they load only through one loader, on the final path and in CI), and the
    scan also reads copilot-validation.md (ids only), the validation's docs and comments, the
    decision records, the change fragments and the changelog. PR #71's skill
    drafts carry the same leak and nothing else. This plan is scanned too, so it points at the
    leak without naming it.
  - **Resolved (2026-10-19, protocol 2026-10-19.2):** the test breathing case was recomposed, not
    the skill: a composition a fresh author chose without dev results, named by no skill (main or
    PR #71's drafts) or snapshot, its seed from a rule fixed before looking. The open-leak list is
    empty; arm A / A+ test scores are unchanged. Approved by the manager on the user's delegation;
    the user may reverse it ([copilot-validation.md](algorithms/copilot-validation.md), Protocol
    2026-10-19.2).

### The A+ ceiling: what the validation can still measure

The engines alone now score 0.99 on dev and on test
(the sealed arm A/A+ test record; [copilot-validation.md](algorithms/copilot-validation.md),
"Batch 1"):

- **Where B can still gain.** The case score lets arm B gain at most 0.01 over A+, on the one
  component A+ misses: the Jahn–Teller long-axis family (`validation/armA.js:171,223` give no
  direction).
- **Why A+ is so high.** The score is dominated by deterministic facts: a perovskite case has
  four symmetry components beside one feature and the spurious-claim test
  (`validation/scoring.js:145-147`), and the point-group facts are exact.

What the live runs can still measure:

1. **B against A** (0.81): the Copilot against the default pages, which is the user-facing
   question.
2. **Whether B matches A+:** reaching the engines' conclusions without dropping one or adding an
   over-claim. A drop below 1 on a case A+ gets right is a model failure the score shows.
3. **What the score does not see:**
   - the gate's send-backs and the evidence still missing at the end;
   - numbers written by reference against typed ones;
   - the computed verdict against the score (earned / over);
   - reproducibility over 3 runs;
   - questions asked;
   - time and tokens.
4. **B against C** on the same measures: what the method adds, the skills together with the
   prompt's method rules, since C removes both (1.8, decided).

Item 1.6 proposes making the comparison with A+ informative on the dev set **(approval)**.

### B and C. Architecture: the deepening candidates and the user's decisions

An architecture review on 2026-10-18 used the improve-codebase-architecture method (Matt Pocock).
It read five areas:

- the agent loop and its guards;
- the tool specs and their dispatch;
- the evidence, reference and verdict modules;
- the engines against their browser ports;
- the harness against the app.

It works in codebase-design's terms: module, interface, depth, seam, adapter, leverage and locality,
plus the deletion test. The plan's earlier sections B (the same science in more than one place) and
C (one tool's facts spread across files) are its evidence. Their findings, such as the duplicated
thresholds, the three shell readings and the side lists, now sit under the candidates below. The
user decided on 2026-10-19:

| # | Deepening | Strength · dependency | Decision |
| --- | --- | --- | --- |
| a | **The case-id leak** (a bug, not a deepening). The harness names the run after the case id (`validation/runFolder.js:30`), which reaches the model as `context.run` in the first message (`src/llm/context/runContext.js:422`). Every arm B and C run is told its case kind. | — | Start, first: **fixed** (protocol 2026-10-19.1) |
| 2 | **The run as the Copilot sees it.** One module turns a run's files and the view state into the Copilot's run. Today App, Dashboard and useAssistant each hold a piece, and the harness copies them (`runFolder.js:15-38`). The browser and the harness become its two adapters. | Strong · local-substitutable | Start, second: **built** (`src/copilotRun.js` with `useCopilotRun.js` and `openRunFolder.js`) |
| 1 | **The turn.** One module runs one question. Today `useCopilot.js:126,137-181` and `harness.js:75-132` each assemble the turn, and they have drifted apart: the harness loses the warmer retry, the no-tools prompt and the tools-on rule, and the gate setup is copied in both. The dock and the harness become its two adapters, and the module emits each request as it was sent. | Strong · mock (model port) | Start, third |
| C | **Arm C's "method off" switch** in the turn module (1.8). | — | Start, after 1: **built** (#98) |
| 3 | **Provider quirks behind the model port.** `streamCopilotChat` takes over the retry temperature, honours or enforces `tool_choice: 'none'`, and reports "tools unsupported" as a typed failure instead of a regex over error text. | Strong · mock | Start, fourth |
| 4 | **One verification module: is this answer verified?** It owns the claim catalogue (phrasing, the analysis that backs a claim, the fact that contradicts it), reference filling, the facts schema, the gate's bookkeeping and the order the checks run in. The prompt's evidence sentence is generated from the catalogue. This replaces the side lists in `evidence.js:86-148`, `verdict.js:61` / `turnState.js:161` and `prompts/copilot.js:32-39`. | Strong · in-process | Start, fifth |
| 7 | **Generated golden fixtures at the orientation and site seams.** The 3σ direction verdict and the split-site statistics are pinned by one generated fixture, not by twin hand-copied pins. | Strong · ports & adapters | Start, in parallel |
| 5 | **The engines return their own verdicts** (clusters, modulations), with the row that decided each, so the tool only renders them. | Strong · ports & adapters | Approved 2026-10-19 (manager, on the user's delegation): A+ results change, and the verdict becomes part of the library output |
| 6 | **One shell reading per pair, per run** (find_gr_peaks, the clarifier, the context, the audit and A+), and a coordination reading in the triplets payload. | Strong · in-process | Approved 2026-10-19 (manager, on the user's delegation): a context flag and the cut rule change |
| 8 | **Arm A+ reads the app's rules** for the split site, the cluster elements and the symmetry floor (A+ uses the floor; the tool and the clarifier use 1.25 × the floor). | Worth exploring · in-process | Approved 2026-10-19 (manager, on the user's delegation): the A+ baseline moves; the app's rule unless there is evidence it is wrong; **built** (#99) |
| 9 | **One request contract for the worker and Flask adapters**: names, defaults, caps, normalisation, error modes, cancelling. | Worth exploring · ports & adapters | Approved 2026-10-19 (manager, on the user's delegation) |
| 10 | **Turn policy: what ends a turn.** One owner for the counters, and one `done` shape. | Worth exploring · in-process | Approved 2026-10-19 (manager, on the user's delegation) |
| 11 | **Tool transcript: what the model is told about a result.** One store for the loop and the UI. Result statuses (old 3.2) attach here. | Speculative · in-process | Approved 2026-10-19 (manager, on the user's delegation) |

**Rejected by the deletion test:**

- **A thresholds constants module (old 2.1).** Its interface would be as wide as its
  implementation: it moves numbers, not decisions. Callers should use the function that already
  decides, for example `isIdealCoordination` (`pairCorrelations.js:130-132`) or
  `facts.direction.significant` (`copilotTools.js:41`).
- **One file per tool in `copilotTools.js` (old 3.4), done before 5 and 6.** It would only spread
  the shallowness over seven files.
- **The two `verdict.js` modules** are different concepts: the validation's case claims and the
  app's computed Outcome. Rename one; don't merge them.

**How each started piece is built:**

- **Design.** Two subagents each propose an interface (design-it-twice) and the lead picks one.
- **PRs.** One PR per piece, stacked a → 2 → 1 → C → 3 → 4. Each targets `main`, so CI runs on it.
  #7 is cut from `main`.
- **What stays fixed.** Each PR keeps the science and the model-visible snapshot (1.7) unchanged,
  with two exceptions: the leak fix (a neutral run name, plus a test that no case id reaches the
  model) and arm C with skills off. Arms A and A+ in `arms.test.js` stay unchanged.
- **Drift between copies.** Where two copies disagree, the PR names the behaviour it keeps and why,
  and pins it with a test.
- **The user.** The user is asked only about choices that change the science, a verdict, a score
  or what the model sees.

**Science items from the old sections B and C that stay open:**

- **The page and the tool read a site at different settings** **(approval, old 2.5).** Arm A reads
  the Displacement Directions page at ν = 10 with smoothing 2 (`armA.js:38-39`); the tool lets the
  engine choose (`copilotTools.js:1149`). On the dev thermal case's O site the map test reads 2.3σ
  against 7.2σ.
- **PCA above 20,000 points draws different subsamples in Flask and in the browser** (old 2.6;
  `rmc_toolkits/pca_kde.py:72,570-576` against `src/workers/pcaKde.js:317-331`). It affects the
  KDE, not a Copilot tool.
- **No test pins the shape of a live tool's output** (old 3.3). With numbers by reference, the data
  paths a result exposes are an interface the model writes against. Closed by 3.3.
- **No corpus of real answers** tests the claim regexes of the checks and the gate (old 3.5).

### D. Tool gaps

These come from the review's tool table that batch 1 left open, plus the survey:

- **`list_sites` leaves out what tells a split site.** The per-axis kurtosis, the axis resolution
  and the PC directions are computed (`pcaKde.js:534`) but not returned (`copilotTools.js:60-70`;
  the review's tool table).
- **No tool returns the static displacement of a split site.**
  - **What the model does now.** It derives δ ≈ √(rms₁² − rms⊥²) itself, a number written by
    hand that provenance flags as unsourced and that numbers by reference cannot carry.
  - **What the scoring accepts.** The scoring's tolerance is built for that estimator
    (`validation/scoring.js:20-21`).
- **The Jahn–Teller bond split is not returned.** `find_gr_peaks` gives the two peaks of a joined
  shell, not their difference (`runTools.js:191`).
- **`analyze_modulations` gives the wavevector, not the displacement direction of the mode.** A+
  takes it from the most prolate site's PC1 (`armA.js:208-211`).
- **Short-range order between cells shows only as a single strong mode.** The general-wavevector
  test (batch 1, L4) reads the strongest modes (`copilotTools.js:836,881`), at most 5 per element
  (`rmc_toolkits/modulations.py:51`). A broad diffuse enhancement over many q is untested.
- **The symmetry tools need a basis** (`copilotTools.js:1303`), which only browser-loaded runs
  have. In Flask mode `detect_symmetry` and `measure_symmetry_lowering` are absent (AGENTS.md,
  next step 4).
- **No validation case has fit data,** so `assess_convergence`, `compare_fit` and the
  rmc-fit-assessment skill are not exercised (the review's skills table).

### E. Skills disagree with the tools

On `main` the skills are the 2026-10-11 drafts: batch 1's skill changes wait in PR #71. Against
the tools on `main`:

- **displacement-anisotropy:**
  - "Preferred direction only with a significance ≥ 3σ" (`SKILL.md:59`) and its list of
    significances (`SKILL.md:34`) count the map test, while the tool calls a direction preferred
    on its peak test only (batch 1, L3).
  - "Name the family only when the angle is small" has no number (`SKILL.md:40`); the shift tool
    uses 15° (`copilotTools.js:978`).
  - The Ta₄ "23 % less distorted" (`significance-reading.md:29`) is 20 % since the shape null.
- **coordination-shells:**
  - The same outdated 23 % (`SKILL.md:85`, `lacunar-spinel-shells.md:19`).
  - Nothing on `first_shell.first_peak_only`, which the tool returns since L1.
- **local-symmetry:**
  - "Which tolerance is the group" is left open (`SKILL.md:57`).
  - `measure_symmetry_lowering` is not mentioned, while the prompt's evidence rule requires it
    for a lowering (`prompts/copilot.js:32-39`).
- **rmc-fit-assessment:** "The verdict is a self-assessment" (`SKILL.md:79`) has been untrue since
  THY-18.
- **No local-order skill.** The method lives in the prompt and the gate instead
  (`prompts/copilot.js:32-39`, `evidence.js:132-148`), a prompt rule that restates a skill.

PR #71 fixes each of these. Three disagreements remain after it:

- **Split sites.** The skill's split site is a negative kurtosis; the tools' is a PC1 kurtosis
  ≤ −0.5 three errors below zero; `list_sites` shows neither (candidates 8 and 4.1).
- **Which elements need clusters.** The local-order draft says "every element whose own first
  shell is short and closed"; the gate says CN ≤ 4.5 and A+ says a peak below 3.6 Å.
- **What "no order" needs.** The local-order draft requires the polyhedra (`find_gr_peaks`) and
  the ladder; the gate's no-order rule requires neither (`evidence.js:132-148`).

Beyond single contradictions:

- **No skill lint.** The skills test checks the sections and the draft line
  (`src/__tests__/copilotSkills.test.js:24-26`). Nothing checks that a named tool exists, that a
  number matches its constant, or that a firm rule the skill claims is in the code.
- **Size against a 32k window.**
  - **Before any analysis** the dev Jahn–Teller question's prompt estimates at about 9.8k tokens
    (the contextGuard's 2.5 characters per token): an 11.0k-character system prompt, 10.0k of
    tool JSON for 11 tools, 2.9k of context.
  - **The other claims on the window.** A skill is 4.7–5.4k characters (about 2k tokens) and is
    kept whole once read (`src/llm/copilot/skills.js:21,121`; `agentLoop.js:618`). Analyses take
    up to 12k characters (about 4.8k tokens). 8,192 tokens are reserved for the reply
    (`contextGuard.js:41`).
  - **The total.** Two skills read then come to about 27k of qwen3.8-27b's 32k: a third skill, or
    a long reasoning round, is refused by the context guard.
  - **Measured (5.5's budget test, #81).** The guard counts tool messages (skills and analyses) at
    1.4 characters per token, not 2.5, so the worst turn does not fit today:
    - **The dock on the Demo:** prompt 4,762 + tool JSON 6,752 + context 2,213 + question 70 + the
      two largest skills 7,814 + analyses at the 12k budget 8,724 = 30,333. With the 8,192
      reserve that is 38,525 tokens, 5,757 over 32,768, or 1,661 over under LM Studio's
      4,096-token rule.
    - **Arm B on its largest dev case:** 26,889 + 8,192 = 35,081, 2,313 over, or 1,783 to spare
      under LM Studio's rule.
    - The numbers grew by 38 tokens with 4.1's description. The fix for the dock is on-demand
      skill references, which is skill text and waits for the user. The live runs load LM Studio
      models at ≥ 48k (Q3).
- **No outline drops the gate numbers.** Unlike MATERIA, the Copilot never shrinks a skill to an
  outline: a read skill stays whole for the question. The risk runs the other way: kept whole, the
  skills crowd the window (above).

## Phases

### Phase 1 — Make the validation measure the Copilot

| PR | Change | Size |
| --- | --- | --- |
| 1.1 | **Done (#84).** **The warmer retry in the harness.** Now part of the turn module (1) and the model port (3): both adapters get the same retry, and a test at the port asserts it. | — |
| 1.2 | **Done (#86).** **One run policy, recorded.** A `RunPolicy` states the limits, temperature, server kind, context size, reply cap and the context-guard decisions; each record and transcript keeps it; the report lists what a comparison across models also compares. One scripted check that the final round's `tool_choice: 'none'` is sent in the form each provider accepts. | M |
| 1.3 | **Done (#78).** **Context parity, stated.** The record says which context sections the run had (no PCA table, no run settings, the 0.2 Å card) and that it is one fresh question; the harness gets an option for a PCA-table context, as after a visit to the PCA page, for a dev check of its effect. **Built with the run module (Phase 2, PR 2):** `record.view`, `record.context`, `run.js live --pca-table`. | S |
| 1.4 | **Done (#90, #93).** **State lifecycle.** A Live Data save that changes the configuration starts a new conversation, as a new run does, or marks earlier answers as about the previous configuration; the Detected SG pick resets with the run. | S |
| 1.5 | **Holdout scan widened** (done, 2026-10-19): test-case materials, formulas, magnitudes and directions (from the test definitions, sealed since #104) across `src`, `.claude/skills`, `rmc_toolkits` and `docs/` except the validation docs. The overlap with the skills' reference compound reaches arm B through `read_skill`: copilot-validation.md lists it as an open leak pending the user. Only test- or doc-only hits are allow-listed, as known risks. | S |
| 1.6 | **An informative comparison with A+** *(approved 2026-10-19; the held-out judgment cases are written by a fresh author without the dev results)*. Dev set only. (i) The report gives feature-only scores, B against A+ per component, and every disagreement with the planted truth beside it, with no change to the case score. (ii) A "judgment" dev subset, designed and tried on dev only: cases where A+'s fixed joining rules are ambiguous or wrong — two features at once, a signal at 2.5–4σ, evidence that two tools read differently, a choice only the user can make — so the model has room to add. The test set gets matching cases only through a new, pinned test version authored without reference to the dev results, decided by the user. | M |
| 1.7 | **Done (#73).** **Snapshots of what the model sees.** A test snapshots the system prompt (every variant), the tool list, every description and schema, the first request as the dock and the harness send it, and a sample of each tool's result text, so a refactor that changes them shows in the diff. Arm C's snapshot is checked for skill names. First of the started set. | S |
| 1.8 | **What arm C ablates: decided 2026-10-19.** C removes the whole method: no skill names, stage lists, `EVIDENCE_RULES`, method sentences in `TOOL_RULES` or method advice in tool descriptions. It is a structural switch in the turn module, after 1 (B and C, item C). | M |

### Phase 2 — The started set (decided 2026-10-19)

The candidates of B and C, in this order. Each is one PR, stacked on the one before.

| PR | Change | Size |
| --- | --- | --- |
| a | **Done (#75).** **The case-id leak.** The harness gives the run a neutral name; a test asserts that no case id or kind reaches the model. The model-visible snapshot changes by the run name only. **(this PR:** `rmc_run`, `caseBlind.test.js`, protocol 2026-10-19.1**)** | S |
| 2 | **Done (#78).** **The run module** (candidate 2). The browser and the harness are its two adapters; the harness records the view state it chose. **(this PR:** `copilotRun.js`, `useCopilotRun.js`, `openRunFolder.js`, the harness on `openRunFolder`, item 1.3**)** | M |
| 1 | **Done (#84).** **The turn module** (candidate 1). The dock and the harness are its two adapters; it emits each request as sent. | L |
| C | **Done (#98).** **Arm C's "method off" switch** in the turn module (1.8). Arm C's snapshot is checked for skill names and the method sentences. | M |
| 3 | **Done (#88).** **The model port** (candidate 3): the retry temperature, `tool_choice: 'none'`, "tools unsupported" as a typed failure. | S |
| 4 | **Done (#101).** **The verification module** (candidate 4). The prompt's evidence sentence is generated from its claim catalogue. | L |
| 7 | **Done (#74).** **Golden fixtures** at the orientation and site seams (candidate 7), in parallel, cut from `main`. Built: `orientation_fixture.json` and `pca_sites_fixture.json` from their generators, with staleness tests on both sides; the hand-copied twin pins are gone. | S |

### Phase 3 — Approved 2026-10-19 (manager, on the user's delegation)

| PR | Change | Size |
| --- | --- | --- |
| 5 | **Done (#89).** **The engines return their own verdicts** (clusters, modulations) *(approved 2026-10-19)*. | M |
| 6 | **Done (#94; the cut rule #113, approved by the user).** **One shell reading per pair, per run** *(approved 2026-10-19)*. Changes the numbers a question quotes (by < 0.05 atom), a context flag and the cut rule. | M |
| 8 | **Done (#99).** **Arm A+ reads the app's rules** for the split site, the cluster elements and the symmetry floor *(approved 2026-10-19)*. It covers old 2.3 and 2.4, which land with 5.3, since the rule and the skill change together. Arm A's rule is part of the frozen protocol and changes only with a dated note. | S |
| 9 | **Done (#92).** One request contract for the worker and Flask adapters: one fixture checked against every adapter. | M |
| 10 | **Done (#103).** Turn policy: what ends a turn. One owner for the counters and one `done` shape. | M |
| 11 | **Done (#107).** Tool transcript; result statuses (old 3.2). | M |
| 2.5 | **Done (#110).** **The page and the tool at the same settings** *(approved 2026-10-19)*: the directions tool also reports at the page's ν and smoothing, or the page defaults to the engine's choice. Changes arm A's reading. | S |
| 2.6 | **Done (#96).** PCA subsamples identical in both runtimes (AGENTS.md next step 8). | M |
| 3.3 | **Done (#77).** **Output-shape contracts** for every live tool: a schema of its `data` (the paths references name), checked on the demo and on each dev case. Before any tool output changes in phase 4. Built: each live tool's spec carries a `dataShape` (`llm/copilot/dataShape.js`, never sent to the model), and `toolDataShapes.test.js` checks every result on the demo and the dev cases against it, both ways (no undeclared key; every declared key returned by some call, the unreachable forms listed). | M |
| 3.5 | A small corpus of real answers for the claim regexes of the checks and the gate. **Built (#83):** `answerCorpus.test.js`, 9 real answers (3 verbatim, 1 transcribed, 5 quoted), one adapter, pins with judgements; 0 false alarms in 4 readings, 1 miss (a negative cluster conclusion is not read). The dev live runs' answers join it before the freeze (6.3). | S |

### Phase 4 — Tool science

Each item ends with a known-answer test on the dev cases or synthetic configurations. Items that
change what A+ sweeps are noted.

| PR | Change | Size |
| --- | --- | --- |
| 4.1 | **`list_sites` shows the split** *(approved 2026-10-19)*: per-axis kurtosis on resolved axes, PC1 as [u v w], and the static displacement of a split site (√(rms₁² − rms⊥²) with its sampling error), so the model quotes it by reference. **Done (#82):** `kurtosis_axes`, `pc1_uvw`, `delta_A` ± `delta_se_A` (ai-assistant.md Step 20); one rule, `isSplitSite`; known answer over 5 × 1,000 seeds (95 % coverage 0.939–0.958; a prolate Gaussian site never called split); arms A/A+ unchanged; protocol 2026-10-19.3. | S |
| 4.2 | **Done (#109).** **The Jahn–Teller split by number** *(approved 2026-10-19)*: `find_gr_peaks` returns the long − short bond of a joined shell and the share in each part. | S |
| 4.3 | **Done (#102).** **The polarisation of an ordered mode** *(approved 2026-10-19)*: `analyze_modulations` returns each ordered point's displacement direction (the mode's eigenvector over the element's sites). Changes A+ (it reads the tool, not PC1). | M |
| 4.4 | **Done (#112).** **Diffuse order between cells** *(approved 2026-10-19; new statistic, calibrated on thermal-only seeds, false-positive rate ≤ α)*: a calibrated test of excess power spread over many q (not one mode), against the same independent-displacement null, with a false-positive test over many seeds of thermal-only configurations. | L |
| 4.5 | **Done (#108).** **The symmetry tools in Flask mode**: the basis (with its moments) from `/api/structure`, so `detect_symmetry` and `measure_symmetry_lowering` exist for server folders (AGENTS.md next step 4). | M |
| 4.6 | **Done (#87).** **A fit-data dev case** for rmc-fit-assessment: a configuration with a synthetic fit and χ² history, so `assess_convergence` and `compare_fit` are exercised. Dev only. | M |

### Phase 5 — Skills

| PR | Change | Size |
| --- | --- | --- |
| 5.1 | **Skill lint in CI.** Every tool a skill names exists; every gate number is bound to the constant of the module that decides it; every "*Firm rule in the app:*" names a check or refusal that exists. It lands with an allow-list of today's violations, which 5.2 to 5.4 empty. **Built (#80):** `src/__tests__/skillLint.test.js` (tables in `skillLint/tables.js`). The allow-list holds three entries, all in displacement-anisotropy: the 15° cone and the peak test (both fixed in #71's drafts) and the split site (5.3). #71 needs a lint update when it merges: 13 bindings, 4 firm-rule mappings and 4 no-gate entries. | M |
| 5.2 | **(approval)** The user's review of PR #71 (batch 1's skill drafts and the new local-order draft). | — |
| 5.3 | **(approval; each draft goes to the user before it merges)** The three disagreements left after #71 (split site, cluster elements, what "no order" needs), in the same PR as candidate 8 and the gate change. | S |
| 5.4 | **(approval; each draft goes to the user before it merges)** Method out of the tool descriptions and the prompt for arm B too: a description says what a tool does; the window advice and the evidence sentence move to the skills. Arm C already drops them structurally (1.8). **Waits for the user's skill review** (manager, 2026-10-19): arm B keeps its method advice until the skill text carries it; arm C drops it through 1.8. | M |
| 5.5 | Size for a 32k window: a budget test (prompt + two skills + the turn's analyses + the reply reserve within 32k); long reference material moves to on-demand references. **The budget test is built (#81):** `contextBudget.test.js`, with today's overflow pinned and the fit an expected failure. The on-demand references wait for the user's skill review. | S |
| 5.6 | SKILL_INTERVIEWS brought up to date when #71 merges; the validation report says which skills are the user's method and which are drafts. | S |

### Phase 6 — Freeze and re-baseline

1. CI green; the scripted replays and the harness tests green.
2. The local tier (qwen3.8-27b, gemma4:26b, meta/muse-glimmer) on the dev set, arms B and C, one
   run per case first (the smoke run), then three: no silent cut, no context refusal, no turn
   lost to the gate.
3. A dry run of the pipeline on the dev set: records, transcripts, scoring, findings, the report's
   THY-18 table.
4. Arms A and A+ again at the freeze; their tables in copilot-validation.md and the pins in
   `arms.test.js` updated.
5. The freeze tagged and its commit written into copilot-validation.md. Then the test set runs
   once, as copilot-validation.md describes. The paid model runs only with the user's go-ahead.

## Order

- **Start together** (2026-10-19): the snapshot test (1.7), then the started set in order (a → 2 →
  1 → C → 3 → 4), with #7 alongside. Phase 1's other items (1.2–1.5) and 5.1 can follow on top.
- **Phase 3** waits for the user's approval batch. 3.3 goes before any tool output changes in phase
  4, so the shape contract shows each change.
- **Phase 4** needs the started set and 3.3.
- **Phase 5's content** (5.3–5.5) waits for phase 4, so the skills describe the final tools; 5.2
  (the review of #71) can happen any time.
- **Phase 6** is last.

The core goes first: 1.7, the started set, 1.2–1.5, 3.3, 4.1 and 5.1–5.3. These change what the
validation would find. The rest is upkeep that follows, before or after the validation.

## Risks

- **A refactor changes behaviour silently.** The snapshot of what the model sees (1.7), the shape
  contracts (3.3), the replays and the engine parity fixtures (with #7's goldens) catch changes in
  what the model sees and in the numbers.
- **The holdout erodes.** Improving the tools for the test cases' kinds is what the dev cases are
  for. The wider scan (1.5) and the dev-only checks keep the test cases themselves unseen; a
  "judgment" test version (1.6) is authored only on the user's decision, without the dev results.
- **The A+ ceiling hides what the model adds.** On the case score alone, B can at best tie A+.
  Unless 1.6 is decided, the report must lead with B against A and the process measures, and say
  so.
- **Local-model robustness regresses.** The loop guards, the gate and the context guard were tuned
  on the local tier. The turn, model-port and verification modules move them without changing
  them; phase 6 runs the local tier again before the test set.
- **The gate shapes the answers it checks.** A send-back costs a round and steers the model toward
  the listed analyses. The record keeps every send-back (THY-18), and 3.5's corpus counts false
  ones before the freeze.
