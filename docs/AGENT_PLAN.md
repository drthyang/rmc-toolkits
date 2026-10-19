# Copilot plan: architecture, tools and skills before the live runs

The validation ([copilot-validation.md](algorithms/copilot-validation.md)) measures the AI Copilot
as it ships. Before arms B and C run on the test set, this plan fixes what would make the numbers
measure the harness instead of the Copilot, puts each piece of science in one place, and brings the
skills into line with the tools. Then the code is frozen, arms A and A+ are run again, and the
validation runs as written.

It builds on the [tool and skill review](algorithms/copilot-review.md) (2026-10-14), validation
batch 1 and THY-18 (PRs #69 and #70, merged), and a survey of `web_app/frontend/src` (the Copilot,
its tools and the workers), `web_app/frontend/validation`, `rmc_toolkits` and `.claude/skills` on
`main` at `b5b3820` (2026-10-18). Every `file:line` below is on that commit. It changes nothing:
the skill drafts of batch 1 wait in PR #71 for the user's review.

## Ground rules

- **The test set stays held out.** Every change here is designed and checked on the dev set, the
  scripted replays and CI. Arms B and C do not run on a test case until the freeze (phase 6). Arms
  A and A+ are deterministic and are re-run at the freeze, since several items change the engines
  A+ sweeps.
- **Approval.** A new tool, a skill change, a tool description or prompt rule that restates a
  skill, and a threshold change that alters the science all need the user's approval. Items waiting
  for it are marked **(approval)**.
- **No API spend.** Local models and scripted CI only, until the user says otherwise. The
  Sonnet arm's $12 share stays unspent.
- **Each PR** passes CI (Python on 3.9 floors, 3.11, 3.13; frontend lint, vitest, build), updates
  its docs, and leaves the app working. A refactor PR keeps the model-visible text unchanged: a
  snapshot of the prompt, the tool list and the descriptions (phase 1.7) shows any change in the
  diff.

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
    0.2 Å (`runFolder.js:15,24-27`).
  - **The prompt and history.** The harness offers no view tools, so the model reads the
    compute-only rules instead of `VIEW_RULES` (`harness.js:90-92`,
    `src/llm/prompts/copilot.js:151`), and it never has earlier turns (`harness.js:91`); the
    app sends up to 8 (`copilot.js:125`, `useCopilot.js:137`). These are known and documented
    choices; the record does not say them.
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
  - **What it checks.** Only the test-case ids, only in `src` (js/jsx), `.claude/skills` and
    `rmc_toolkits` (`validation/__tests__/split.test.js:53-59`).
  - **What it misses.** The test cases' materials, formulas, magnitudes and directions, `docs/`,
    and the llm README.
  - **Today there are no hits;** nothing keeps it that way. One overlap is by design and worth
    stating: the test breathing case's motif (Nb₄ tetrahedra in a selenide,
    `validation/cases.js:509`) is the skills' reference compound, whose published Nb₄ distortion
    they describe (`.claude/skills/coordination-shells/lacunar-spinel-shells.md:31-32`). A model
    may expect a distortion there.

### The A+ ceiling: what the validation can still measure

The engines alone now score 0.99 on dev and on test
(`validation/baseline/arm-a-test.json`; [copilot-validation.md](algorithms/copilot-validation.md),
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
4. **B against C** on the same measures (what the skills add on top of a prompt that already
   carries the method, above).

Item 1.6 proposes making the comparison with A+ informative on the dev set **(approval)**.

### B. The same science in more than one place

The engines exist twice by design: `rmc_toolkits/*.py` is the source of truth, and the browser
ports in `web_app/frontend/src/workers/` are parity-tested against Python goldens. The survey found:

- **PCA above 20,000 points differs between runtimes.** Python draws the fit subsample with numpy's
  generator (`rmc_toolkits/pca_kde.py:72,570-576`), the port with mulberry32 and Fisher–Yates
  (`src/workers/pcaKde.js:317-331`, which says so).
  - **Effect.** A pooled element cloud, or a site in a box of 28 cells per edge or more, gives
    different numbers in Flask mode and in the browser.
  - **Which runtime the validation measures.** Only the browser runtime (`runFolder.js:28`); the
    Flask path, which serves typed server folders, is untested by it.
- **The page and the tool read the same site differently.**
  - **The settings.** Arm A reads the Displacement Directions page at ν = 10 with smoothing 2
    (`armA.js:38-39`). The tool lets the engine choose the resolution (`copilotTools.js:1149`).
  - **The effect.** On the dev thermal case's O site the map test reads 2.3σ at ν = 10 and 7.2σ at
    the engine's choice, so "anisotropic" holds for the tool and not for the page.
  - **The skill** notes the difference ("say which you quote") without a rule.
- **Three readings of one g(r) shell.**
  - **The readings.** `find_gr_peaks` (`src/llm/copilot/runTools.js:179,191`), the Workbench's
    shell catalogue (`src/llm/copilot/clarify.js:91-120`) and the context's `pair_correlations`
    (`src/llm/context/pairCorrelations.js`) all call `firstShell`, then round and integrate apart:
    - the window: 3 significant figures against 0.01 Å;
    - the CN: grid trapezoid to the boundary index (`pairCorrelations.js:251`) against
      interpolated trapezoid over the rounded window (`clarify.js:54,91`);
    - the "first peak" window: r ± FWHM (`runTools.js:179`) against the joined shell's own
      boundary (`clarify.js:119`).
  - **The effect.** A question's option can quote a CN that differs from the one `find_gr_peaks`
    shows the model.
- **Three definitions of a split site:**
  - Arm A: prolate ≥ 1.3, PC1 resolved and PC1 excess kurtosis ≤ −0.5 (`armA.js:35-36,77`);
  - the cluster summary: PC1 kurtosis ≤ −0.5 and 3 sampling errors below 0
    (`copilotTools.js:489-497`);
  - `list_sites`, which the model reads: Mardia's kurtosis only (`copilotTools.js:68`), although
    the engine computes the per-axis kurtosis and resolution (`src/workers/pcaKde.js:534`).

  So the model cannot apply the rule its own tools apply.
- **Two definitions of an element that forms clusters.** The A+ sweep tries clusters for a
  same-element first shell peaking below 3.6 Å (`armA.js:179`). The evidence gate requires
  `analyze_clusters` for a first shell of at most 4.5 neighbours (`evidence.js:61`). A case can
  need clusters by one rule and not the other.
- **Thresholds defined more than once:**

  | Quantity | Copies |
  | --- | --- |
  | The Detected SG card's 0.2 Å | `App.jsx:92`, `validation/runFolder.js:15`, `armA.js:47`, `copilotTools.js:1313,1341,1344` |
  | 1.25 × the noise floor | `clarify.js:50` (`NOISE_MARGIN`), `copilotTools.js:29` (`NOISE_FLOOR_MARGIN`) |
  | The long-wave limit 0.25 | `rmc_toolkits/modulations.py:58`, `src/workers/modulations.js:23`, `copilotTools.js:779` |
  | A shell not separated: g at 2 % of the peak | `pairCorrelations.js:273`, `runTools.js:137`, `clarify.js:178` (a literal) |
  | Half an atom, the CN tolerance | `pairCorrelations.js:130`, `clarify.js:279,289`, `copilotTools.js:322,329`, `src/llm/copilot/answerAudit.js:378` (literals) |
  | 3σ | `copilotTools.js:41`, `armA.js:40`, `prompts/copilot.js:54` and `answerAudit.js:241` (in text) |
  | Split-site kurtosis −0.5 | `armA.js:36`, `copilotTools.js:489` |
  | 15° | `rmc_toolkits/clusters.py:81` and `src/workers/clusters.js:29` (axis cone), `copilotTools.js:978` (low-index cone), `validation/scoring.js:28` (scoring) |

  Only the engine pairs are pinned equal by a fixture. The rest can drift apart unseen.

### C. One tool's facts are spread across files

- **Per-tool facts sit in side lists.** The tool spec carries its page, cards, page parameters,
  bond windows and result budget (`copilotTools.js`). Several other facts live elsewhere:

  | List | File |
  | --- | --- |
  | which tool evidences which conclusion | `evidence.js:86-148` and, again in prose, `prompts/copilot.js:32-39` |
  | which claim kind maps to which answer check | `validation/scoring.js:154` |
  | which steps are analyses | `verdict.js:61` (no views, no unrun calls) against `turnState.js:161` (counts views) |
  | the window question's tool | `harness.js:234` (by name) |

  The analysis count already disagrees: the footer counts a view tool as an analysis
  (`src/llm/copilot/turnState.js:161`), the computed verdict does not (`src/llm/copilot/verdict.js:61`).
- **Results are read back as text** in two places:
  - an unknown tool, by matching its error message (`harness.js:240`);
  - a kept window, by the option's label starting with "Keep" (`clarify.js:510`).

  The answer checks and the gate read prose by regex by design. Those regexes have no corpus test
  over real answers beyond the replays.
- **One effect label is loose.** `orient_3d_view` is a view tool but selects a site as a click
  would (`src/copilotViewTools.js:9,60`), which sets the input of the PCA page's own analysis. That
  was the maintainer's decision (AGENTS.md, view tools); it should be visible in the spec.
- **No test pins the output shape of the live tools.** With numbers by reference, a result's data
  paths are an interface the model writes against (`{{r3.sites.0.shift_A}}`). The tests pin chosen
  fields (`src/__tests__/copilotTools.test.js`), not the shape, so a renamed key silently breaks
  every reference to it.
- **The biggest files mix concerns:**
  - `copilotTools.js` (1,406 lines) holds eight summarizers, the tool specs and the React hook;
  - `agentLoop.js` (650) holds the loop, the clarifier calls, compaction, the gate, the reference
    fill and the verdict;
  - `answerAudit.js` (430) holds claim reading and ten rules;
  - `clarify.js` (514) holds the shell catalogue and the questions.

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
  ≤ −0.5 three errors below zero; `list_sites` shows neither (B above).
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
- **No outline drops the gate numbers.** Unlike MATERIA, the Copilot never shrinks a skill to an
  outline: a read skill stays whole for the question. The risk runs the other way: kept whole, the
  skills crowd the window (above).

## Phases

### Phase 1 — Make the validation measure the Copilot

| PR | Change | Size |
| --- | --- | --- |
| 1.1 | **The warmer retry in the harness.** `run.js` honours `antiRepeat` as `useCopilot` does (0.6 at least). A harness test asserts it. | S |
| 1.2 | **One run policy, recorded.** A `RunPolicy` states the limits, temperature, server kind, context size, reply cap and the context-guard decisions; each record and transcript keeps it; the report lists what a comparison across models also compares. One scripted check that the final round's `tool_choice: 'none'` is sent in the form each provider accepts. | M |
| 1.3 | **Context parity, stated.** The record says which context sections the run had (no PCA table, no run settings, the 0.2 Å card) and that it is one fresh question; the harness gets an option for a PCA-table context, as after a visit to the PCA page, for a dev check of its effect. | S |
| 1.4 | **State lifecycle.** A Live Data save that changes the configuration starts a new conversation, as a new run does, or marks earlier answers as about the previous configuration; the Detected SG pick resets with the run. | S |
| 1.5 | **Holdout scan widened**: test-case materials, formulas, magnitudes and directions (from `cases.js`) across `src`, `.claude/skills`, `rmc_toolkits` and `docs/` except the validation docs. The overlap with the skills' reference compound is listed in copilot-validation.md as a known risk. | S |
| 1.6 | **An informative comparison with A+** **(approval)**. Dev set only. (i) The report gives feature-only scores, B against A+ per component, and every disagreement with the planted truth beside it, with no change to the case score. (ii) A "judgment" dev subset, designed and tried on dev only: cases where A+'s fixed joining rules are ambiguous or wrong — two features at once, a signal at 2.5–4σ, evidence that two tools read differently, a choice only the user can make — so the model has room to add. The test set gets matching cases only through a new, pinned test version authored without reference to the dev results, decided by the user. | M |
| 1.7 | **Snapshots of what the model sees.** A test snapshots the system prompt (both arms), the tool list, every description and schema, and the gate's send-back text, so a refactor that changes them shows in the diff. Arm C's snapshot is checked for skill names. | S |
| 1.8 | **What arm C ablates, decided** **(approval)**. Either C stays "B without the skills' text" and the doc says so, or C also drops `EVIDENCE_RULES`, the method sentences of `TOOL_RULES` and the tool descriptions' window advice. That makes C "no method" and needs a structural switch, not text edits. | M |

### Phase 2 — One definition per rule

| PR | Change | Size |
| --- | --- | --- |
| 2.1 | **A thresholds module** (`src/thresholds.js`): the card tolerance, the noise-floor margin, the long-wave limit (with the Python constant pinned by the fixture), the separation fraction, the CN tolerance, the significance level, the split-site kurtosis and the cones, imported by the tools, the clarifier, the checks, the gate and arm A. No value changes; the skill lint (5.1) reads it. | M |
| 2.2 | **One shell reading.** `find_gr_peaks`, the clarifier's catalogue and the context take the window, CN and first-peak window from one function, with one rounding. Changes the numbers a question quotes (by < 0.05 atom). | M |
| 2.3 | **One split-site rule** **(approval)**: the per-axis rule (PC1 resolved, excess kurtosis ≤ −0.5 and three sampling errors below 0) used by the cluster summary, arm A+ and, through 4.1, shown to the model. Arm A's rule is part of the frozen protocol and changes only with a dated note. | S |
| 2.4 | **One rule for which elements need clusters** **(approval)**: the gate's CN ≤ 4.5, A+'s 3.6 Å and the local-order draft's wording made one. It changes what the gate requires and what A+ sweeps. | S |
| 2.5 | **The page and the tool at the same settings** **(approval)**: the directions tool also reports at the page's ν and smoothing, or the page defaults to the engine's choice (ROADMAP: Displacement Directions default resolution). Changes arm A's reading. | S |
| 2.6 | **PCA subsamples identical in both runtimes**: the Python draw ported to the browser or the other way (AGENTS.md next step 8). | M |

### Phase 3 — One tool spec, structured results

| PR | Change | Size |
| --- | --- | --- |
| 3.1 | The tool spec carries every per-tool fact: what it evidences (`evidence.js` reads it; the prompt's evidence sentence is generated from it), whether it is an analysis, view or question, and the claim kinds its facts cover (`scoring.js` reads it). A test fails if a consumer keeps its own list. The analysis count is one helper. | M |
| 3.2 | Results carry a status (ok, refused, error, unknown tool, skipped) and a choice flag on clarifications; the harness and the clarifier read them instead of text. | S |
| 3.3 | **Output-shape contracts** for every live tool: a schema of its `data` (the paths references name), checked on the demo and on each dev case. | M |
| 3.4 | `copilotTools.js` splits into one module per tool (summarizer + spec) and the hook; `agentLoop.js` into the loop, the clarifier step and the finish (references, checks, gate, verdict). No behaviour change; the replays are the net. | L |
| 3.5 | A small corpus of real answers (the 2026-10-11 spot check, the dev transcripts once they exist) for the claim regexes of the checks and the gate: each finding and each send-back pinned, false alarms counted. | S |

### Phase 4 — Tool science

Each item ends with a known-answer test on the dev cases or synthetic configurations. Items that
change what A+ sweeps are noted.

| PR | Change | Size |
| --- | --- | --- |
| 4.1 | **`list_sites` shows the split** **(approval)**: per-axis kurtosis on resolved axes, PC1 as [u v w], and the static displacement of a split site (√(rms₁² − rms⊥²) with its sampling error), so the model quotes it by reference. | S |
| 4.2 | **The Jahn–Teller split by number** **(approval)**: `find_gr_peaks` returns the long − short bond of a joined shell and the share in each part. | S |
| 4.3 | **The polarisation of an ordered mode** **(approval)**: `analyze_modulations` returns each ordered point's displacement direction (the mode's eigenvector over the element's sites). Changes A+ (it reads the tool, not PC1). | M |
| 4.4 | **Diffuse order between cells** **(approval; new statistic)**: a calibrated test of excess power spread over many q (not one mode), against the same independent-displacement null, with a false-positive test over many seeds of thermal-only configurations. | L |
| 4.5 | **The symmetry tools in Flask mode**: the basis (with its moments) from `/api/structure`, so `detect_symmetry` and `measure_symmetry_lowering` exist for server folders (AGENTS.md next step 4). | M |
| 4.6 | **A fit-data dev case** for rmc-fit-assessment: a configuration with a synthetic fit and χ² history, so `assess_convergence` and `compare_fit` are exercised. Dev only. | M |

### Phase 5 — Skills

| PR | Change | Size |
| --- | --- | --- |
| 5.1 | **Skill lint in CI.** Every tool a skill names exists; every gate number is bound to its constant in the thresholds module (2.1); every "*Firm rule in the app:*" names a check or refusal that exists. It lands with an allow-list of today's violations, which 5.2 to 5.4 empty. | M |
| 5.2 | **(approval)** The user's review of PR #71 (batch 1's skill drafts and the new local-order draft). | — |
| 5.3 | **(approval; each draft goes to the user before it merges)** The three disagreements left after #71 (split site, cluster elements, what "no order" needs), in the same PR as 2.3, 2.4 and the gate change. | S |
| 5.4 | **(approval; each draft goes to the user before it merges)** Method out of the tool descriptions and the prompt: a description says what a tool does; the window advice and the evidence sentence move to the skills, together with 1.8. | M |
| 5.5 | **(approval)** Size for a 32k window: a budget test (prompt + two skills + the turn's analyses + the reply reserve within 32k); long reference material moves to on-demand references. | S |
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

- **Start together.** Phase 1 (1.1–1.5, 1.7) and 5.1 start together. 1.6 and 1.8 wait for the
  user's decisions.
- **Then phase 2,** with 2.1 first. 2.3 and 2.4 land with 5.3, since the rule and the skill change
  together.
- **Phase 3** can overlap phase 2 once 2.1 is in; 3.3 goes before any tool output changes in
  phase 4, so the shape contract shows each change.
- **Phase 4** needs phases 2 and 3.3.
- **Phase 5's content** (5.3–5.5) waits for phase 4, so the skills describe the final tools; 5.2
  (the review of #71) can happen any time.
- **Phase 6** is last.

The core goes first: 1.1–1.5, 1.7, 2.1–2.4, 3.3, 4.1, 5.1–5.3. These change what the validation
would find. The rest is upkeep that follows, before or after the validation.

## Risks

- **A refactor changes behaviour silently.** The snapshot of the model-visible text (1.7), the
  shape contracts (3.3), the replays and the engine parity fixtures catch changes in what the model
  sees and in the numbers.
- **The holdout erodes.** Improving the tools for the test cases' kinds is what the dev cases are
  for. The wider scan (1.5) and the dev-only checks keep the test cases themselves unseen; a
  "judgment" test version (1.6) is authored only on the user's decision, without the dev results.
- **The A+ ceiling hides what the model adds.** On the case score alone, B can at best tie A+.
  Unless 1.6 is decided, the report must lead with B against A and the process measures, and say
  so.
- **Local-model robustness regresses.** The loop guards, the gate and the context guard were tuned
  on the local tier. 3.4 moves them without changing them; phase 6 runs the local tier again
  before the test set.
- **The gate shapes the answers it checks.** A send-back costs a round and steers the model toward
  the listed analyses. The record keeps every send-back (THY-18), and 3.5's corpus counts false
  ones before the freeze.
