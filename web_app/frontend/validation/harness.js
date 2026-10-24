// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// The live arms: the Copilot itself on a case, minus the screen. Each run is
// one question through the turn module (src/llm/copilot/turn.js), the same
// assembly the dock runs (CopilotChat → useCopilot): the context, the prompt,
// ask_user and read_skill, the Workbench's questions, the evidence gate, the
// answer checks and the model port's warmer retry. This file is that module's
// second adapter: it adds the case's run, its tools, the simulated user, the
// arm's skills and the port, and keeps the record.
//
//   B — tools, skills (read_skill + the skills' index), the method in the
//       prompt, the tool descriptions and the checks, guardrails (tool
//       refusals, the Workbench's window questions, the evidence gate, the
//       answer checks);
//   C — the same without the method: the ablation. The turn's method 'off'
//       (src/llm/copilot/method.js, plan item 1.8): no skills, no
//       EVIDENCE_RULES, no method sentence in the prompt's rules, no method
//       advice in a tool description or a check; the guardrails still run.
//
// The view tools (highlight_cards, orient_3d_view) need a screen and are left
// out of both arms; every tool offered only computes. The simulated user
// skips every question (the model must use its stated default), and the run
// and its files have one neutral name for every case (runFolder.js RUN_NAME:
// a case id names its kind), so nothing about the planted truth reaches the
// model except through the tools. The record keeps the case id; the model
// never sees the record.
//
// The run is opened through the run module's in-process adapter
// (runFolder.js loadRun → src/openRunFolder.js), the app's own path from a
// folder's files to the Copilot's run props, tools and engine, with a view
// (HARNESS_VIEW: the app's defaults). The record keeps the view and the
// context sections the model was given, and that the question was one fresh
// question with no earlier turns; `pcaTable` adds the PCA table, as after a
// visit to the PCA Ellipsoid page, for a dev check of its effect (off by
// default).
//
// One run gives one record: every model request as it was sent, every
// decision (round, tool, arguments, the reasoning and text the model wrote
// before it, the result it read), how the turn ended (`stopped_by`, from its
// terminal event: a key of the turn policy's STOPPED_BY, as the dock reads
// it), the final answer, its verdict and score,
// the answer checks against the contradictions the score found, where each
// reported number came from (against what the model was told), and the
// safety counts. The record is the transcript's source.
//
// Each step is its tool_end as the loop sent it (the tool transcript,
// architecture candidate #11): beside the severity (`status`) it keeps how
// the call ended (`result_status`, the tool contract's: ok, refused, error,
// unknown_tool or skipped; src/llm/copilot/tools.js), the note the app added
// for the model (`note`), the user's reply to a question (`choice`), the
// result's reference (`ref`) and whether it is an analysis (`analysis`). A
// step still running when the turn was stopped or failed is listed too,
// `status: 'stopped'` with `result_status` null; no count reads it. The
// safety counts are read from `result_status`, never from the error's text.
//
// THY-18 (protocol 2026-10-17.2): the evidence gate is on (an answer whose
// analysis did not run goes back to the model; each send-back is recorded),
// numbers may be written by reference (the record keeps the answer as written
// and as filled; the VERDICT is parsed from the filled text), and the Outcome
// is the app's computed verdict, recorded beside the case score. All three
// are the turn's verification (src/llm/copilot/verification/); the harness
// adds the verdict's words to the claims it reads (verdict.js
// VERDICT_PHRASINGS, the turn's `config.phrasings`), which the app's own
// catalogue does not read.
//
// The run policy (plan item 1.2, protocol 2026-10-19.5; runPolicy.js): the
// record keeps the rules the run ran under (`record.policy`: the server and
// the window it read before the run, the limits, the temperature, the reply
// cap, the context guard's decisions, the final round's tool_choice form),
// and each request the num_ctx Ollama was sent. validation/run.js runs the
// pre-run check (runPolicy.js preflight) and passes what it read as `server`.
//
// The model port (protocol 2026-10-19.7; src/llm/provider/model.js): each
// request in `requests` is the request as sent with the context guard's plan
// (`sent.plan`: the server kind, Ollama's numCtx / numPredict, LM Studio's
// loaded context, the tool choice's form); a request the guard refuses is
// never sent, has no entry, and the run's `error` is the ContextLimitError.

import { ASK_USER } from '../src/llm/copilot/askUser';
import { READ_SKILL } from '../src/llm/copilot/skills';
import { buildRunTools } from '../src/llm/copilot/runTools';
import { ERROR_STATUSES } from '../src/llm/copilot/tools';
import { runTurn } from '../src/llm/copilot/turn';
import { HARNESS_VIEW, loadRun } from './runFolder';
import { parseVerdict, VERDICT_INSTRUCTIONS, VERDICT_MARK, VERDICT_PHRASINGS } from './verdict';
import { outcomeAgreement, scoreVerdict } from './scoring';
import { shownValues, traceNumbers } from './provenance';
import { buildRunPolicy, noteContextUsed } from './runPolicy';

// Bumped when the cases, question, scoring, harness, tools or skills change what
// a run means (2026-10-17.1: validation batch 1's tools; 2026-10-17.2:
// THY-18's evidence gate, numbers by reference and computed verdict;
// 2026-10-19.1: the run is named RUN_NAME, no longer after the case id, which
// told the model the case's kind; 2026-10-19.2: the test breathing case recomposed,
// its old cluster being the skills' reference compound's, which read_skill
// serves to arm B; 2026-10-19.3: list_sites shows the split — per-axis kurtosis,
// PC1 as [u v w], a split site's static displacement — and says so in its
// description;
// 2026-10-19.4: runs through the turn module, as the dock:
// the request asked again after a loop at a temperature of 0.6 or more, the
// no-tools prompt for a model that refuses tools, provenance against the
// context as shown, every request recorded;
// 2026-10-19.5: the run policy — an LM Studio model runs only when loaded
// with at least 48k, checked before every run; the record keeps the policy
// and the window it ran under; what the model is sent is unchanged;
// 2026-10-19.6: dev-fit-unconverged joins the dev set, the first case with fit
// files (an F(Q), its G(r) and a χ² log; plan item 4.6);
// 2026-10-19.7: the model port — an Ollama run's final request carries no
// tools (tool_choice 'none' enforced; the native API has none), each request
// as sent keeps the context guard's plan, and a request the guard refuses is
// not recorded as sent;
// 2026-10-19.8: one shell reading per pair, per run (src/llm/context/
// shells.js) — the first peak alone ends at the minimum between the first
// two resolved peaks, no longer r ± FWHM; the context's pair_correlations and
// find_gr_peaks read peaks below 8 Å and integrate each CN over exactly the
// window they quote; find_gr_peaks' first_peak_window_A leaves for
// first_shell.first_peak_only, and its shell no longer depends on r_max;
// 2026-10-19.9: the judgment subset joins the dev set (plan item 1.6 ii:
// cases where A+'s fixed joining rules are ambiguous or wrong; dev only until
// a fresh author writes their held-out matches), so a dev batch now holds
// them, and its means are not an earlier batch's;
// 2026-10-19.10: arm C follows no method (the turn's method 'off', plan item
// 1.8) — no EVIDENCE_RULES, no method sentence in TOOL_RULES or ASK_RULES, the
// domain rules' definitions and conduct without their inference, no method
// advice in a tool description or a check, and an evidence-gate send-back
// that names no analysis; arm B unchanged;
// 2026-10-19.11: five held-out judgment cases join the test set, written
// blind from judgment-spec.md;
// 2026-10-19.12: analyze_modulations gives each ordered point's polarisation
// (plan item 4.3), never prints a ratio or σ below its bar rounded onto it
// (F2), and corrects over every element's points whichever the call names
// (F3); arm A+ reads a modulation's direction from the polarisation);
// 2026-10-19.13: the tool transcript's stricter safety reading (architecture
// candidate #11, drift D3) — "asked before" counts only a question that
// reached the user (an ask_user step with the user's reply as its choice),
// not a malformed or refused one, nor a window the Workbench settled with
// no card, so unsafe_actions can rise where a call marked as the user's
// followed only such a question; and (drift D6) number provenance credits
// only what the model was told — each tool message it read, its error and
// note included, never a Workbench question's card — so a number only the
// card held is unsourced, and one only a note or an error held is sourced;
// what the model is sent is unchanged;
// 2026-10-19.15: find_gr_peaks gives a shell joined across a deep trough its
// split by number (plan item 4.2) — each part's coordination number, share
// and mean distance, and long − short with its sampling error, in
// first_shell.split, its description and its joined-shell note (which no
// longer quotes the peaks' separation); arm A+ reads a Jahn–Teller claim's
// magnitude from it;
// 2026-10-19.16: the modulation engines test against independent cells, each
// cell's own correlations kept (F1: sites sharing a cell-level variable no
// longer read as order, and an independent site is no longer held to the
// one-component bound, so both dev weak orders are above the bar), and give
// each element's diffuse order between cells (plan item 4.4), which
// analyze_modulations says and A+ does not read;
// 2026-10-19.17: the Displacement Directions page and the tool at the same
// settings (plan item 2.5) — the page opens at the engine's resolution (Auto)
// and reads every readout from the raw histogram, so arm A (and A+, through
// the same reading) reads a split site's direction map at the tiling
// analyze_displacement_directions reads, from the page's own defaults
// (src/directionSettings.js); and the direction resolution, a measurement
// limit — the engines give the peak cell's circumscribed radius, the tool
// gives the peak direction ± that radius and, above 15° (LOW_INDEX_CONE_DEG),
// says it is not resolved finer and names no low-index family (its
// description and result change; nothing else the model is sent does), arm
// A+ takes PC1 over a significant peak whose cell is wider than 15°
// (directionResolved), arm A reads the page as a user would, and both record
// the radius; the page's map opens unsmoothed (UI only);
// 2026-10-19.18: compute_bond_angles refuses on the one cut rule, the
// Workbench's (shells.js windowCut, approved by the user): a window that meets
// its first shell is cut when it holds fewer partner atoms than the shell by
// at least max(0.5, 15 % of its CN), whatever its coordination counts, which
// decide only where the run reads no shell to compare (as before); the
// refusal names the CN shortfall by path (shell_ab.cn_in_window against
// shell_ab.cn), and every result adds the shell reading (shell_ab, shell_bc)
// and the counts' modal, modal_share and unbonded_share. A whole shell is
// never refused (dev-offcentre's Nb–Nb), a sub-shell first peak always unless
// the user chose it (the Demo's Ta–Se, dev-jahn-teller's Cu–F);
// 2026-10-19.19: the marginal situation re-set for the independent-cell null
// — the dev weak-order pair's amplitude set so its expected (median)
// significance at q is 3σ (judgment-spec.md's new formula), its seeds by its
// rule, one on each side of the bar again — and arm C reads analyze_modulations'
// long-wave and in-phase notes without what the pattern is evidence of
// (inference); arm B unchanged.
// 2026-10-19.20: the held-out marginal pair re-set by a fresh blind author
// from the updated spec (the marginal situation was re-set after F1's
// corrected null raised the test's power, 2026-10-19).
// 2026-10-24.1: the user's fit-assessment interview (approved by the user,
// 2026-10-24) — convergence judged on every χ² term of the .log and their
// total, converged only when none still falls beyond 3 × its noise over the
// last tenth (assess_convergence's description and result, the context's
// convergence.last_10pct in place of change_over_last_10pct and
// recent_slope_per_step, the system prompt's watchdog_status sentence; the
// context budget 6,500 characters); find_gr_peaks calls a peak
// constraint-limited above 5 % of the peak at the cut, and advises (method)
// interpreting nothing on its low-r side and rerunning with the constraint
// loosened; and two checks on the answer, order claimed at a box-scale
// wavevector (long-wave) and a distribution's numbers quoted without saying
// they come from one configuration (single-configuration). Arms A and A+ on
// dev are unchanged.
// 2026-10-24.2: the user's interview answers for four skills (approved by the
// user, 2026-10-24), as checks on the answer and on the tools' results —
// coordination-shells: every CN or bond angle in an answer carries its window
// in Å (window), the coordination check reads "n A per B" as B's count and a
// clause stating the whole shell's window as a claim about the shell, and
// find_gr_peaks cross-checks first_shell.window_A against the average
// structure's bond lengths for the pair (first_shell.average_structure_d_A;
// a warning past the first peak's FWHM; its description); the evidence gate
// backs distorted polyhedra with polyhedral_distortion where the run offers
// it, else analyze_clusters (EVIDENCE_RULES; an anyOf row asks for the first
// analysis the turn offers), and analyze_clusters' method keeps it to finite
// clusters; local-symmetry: a group–subgroup relation that does not hold
// (subgroup: against the ladder's nested groups — index, t or k — or the
// point groups alone) and a cell doubling the box cannot hold (cell-doubling:
// odd along an axis it doubles, or along axes no ordered zone-boundary point
// doubles); detect_symmetry's facts give its operations and tolerance; the
// CIF calls itself an average-structure model of an RMC configuration.
export const PROTOCOL_VERSION = '2026-10-24.2';

// `method` is the turn's (src/llm/copilot/method.js); an arm without one
// follows the full method.
export const ARMS = Object.freeze({
    B: { skills: true, label: 'Copilot with skills and guardrails' },
    C: { skills: false, method: 'off', label: 'Copilot without the method: no skills, method rules or advice (ablation)' }
});

/** The simulated user: skips every question, so the model uses its own default. */
export const skippingUser = (log) => async (request) => {
    log.push({ question: request.question, options: request.options.map((option) => option.label), reply: 'skipped' });
    return { source: 'skipped' };
};

const clip = (text, limit) => (text && text.length > limit ? `${text.slice(0, limit)}… [${text.length - limit} more characters]` : text);

/** The question each live run is asked: the case question plus the verdict format. */
export const liveQuestion = (built) => `${built.question}\n\n${VERDICT_INSTRUCTIONS}`;

// The sections buildRunContext (src/llm/context/runContext.js) can give the
// model besides the run's name and Live Data flag, in its order. A run in the
// app can have them all (the Demo with the PCA table does: runFolder.test.js);
// a case has the structure, its symmetry and the partial g(r), unless the view
// adds the PCA table; a case with fit files (runFolder.js caseFiles) also has
// the optimisation, the datasets, the partial weights and the convergence.
export const CONTEXT_SECTIONS = Object.freeze([
    'structure', 'configuration_optimization', 'run_settings', 'symmetry', 'pca_displacements',
    'pair_correlations', 'datasets', 'partial_weights', 'convergence'
]);

/**
 * What the model was told about the run, for the record: the context sections
 * it had and lacked, and that the question came with no earlier turns (the
 * app sends up to eight).
 */
export const contextRecord = (context) => ({
    sections: CONTEXT_SECTIONS.filter((key) => key in context),
    absent: CONTEXT_SECTIONS.filter((key) => !(key in context)),
    other_sections: Object.keys(context).filter((key) => key !== 'run' && key !== 'live_mode' && !CONTEXT_SECTIONS.includes(key)),
    earlier_turns: 0
});

/**
 * One live run. `port` is the model port the run is sent through
 * (src/llm/provider/model.js: validation/run.js gives providerModel bound to
 * the model, metered for an API model; tests give a scripted fake); `skills`
 * is `{ files, references, load }` as copilotSkills.js gives the app (read in
 * arm B only). `view` (merged over HARNESS_VIEW) and `pcaTable` set what is
 * "on screen" (runFolder.js loadRun); the defaults change nothing. `server`
 * is what the pre-run check read (runPolicy.js preflight; null: not read): it
 * goes into the record's run policy, with the num_ctx each request was sent
 * (read from the context guard's plan in the request as sent).
 */
export const runLiveCase = async ({
    built, arm, model, runIndex = 0, port, skills = null, commit = null, limits = {}, signal = null, now = () => Date.now(),
    view = HARNESS_VIEW, pcaTable = false, server = null
}) => {
    if (!ARMS[arm]) throw new Error(`unknown arm '${arm}' (B or C)`);
    const startedAt = now();
    const run = await loadRun(built, view, { pcaTable });
    const question = liveQuestion(built);
    const asked = [];

    const record = {
        id: `${built.id}__${arm}__${model.name}__r${runIndex + 1}`,
        protocol: PROTOCOL_VERSION,
        commit,
        case: built.id,
        split: built.split,
        arm,
        arm_label: ARMS[arm].label,
        model: { name: model.name, label: model.label, id: model.model, base_url: model.baseUrl, temperature: model.temperature },
        run_index: runIndex + 1,
        started_at: new Date(startedAt).toISOString(),
        // The rules the run ran under (runPolicy.js): what a comparison
        // across models also compares.
        policy: buildRunPolicy({ model, limits, server }),
        question,
        // What the turn offered (its `start`): the method it followed, the
        // tools and the skills.
        method: null,
        tools_offered: [],
        skills_offered: [],
        // The view the run was opened with (the run module's, normalised) and
        // what the model was told: its context sections, one fresh question.
        view: run.view,
        context: null,
        // Every model request as it was sent (the turn's `request` events):
        // { round, why, sent }, `sent` as the port prepared it.
        requests: [],
        steps: [],
        retracted: [],
        gates: [],
        rounds: 0,
        final: null,
        error: null,
        // How the turn ended, from its terminal event (src/llm/copilot/
        // turnPolicy.js STOPPED_BY): answer, max_rounds, no_progress, empty,
        // stopped or error.
        stopped_by: null
    };

    // The round's reasoning and text, the "why" of the calls that follow it.
    let round = { index: 0, reasoning: '', text: '' };
    const open = new Map();
    let gated = null;
    // The run context as the model read it (the turn's `start`).
    let shown = null;
    // One question through the turn module (src/llm/copilot/turn.js), as the
    // dock runs it; the harness adds only what is its own: the case's run, its
    // compute-only tools (the tools over the run props, with no cards, and the
    // run's engine tools), the simulated user, the arm's skills and method,
    // the verdict's words for the evidence gate, and its port.
    for await (const event of runTurn({
        question,
        run: run.props,
        model: port,
        tools: [...buildRunTools(run.props, { cards: [] }), ...run.tools],
        askUser: skippingUser(asked),
        skills: ARMS[arm].skills ? skills : null,
        signal,
        clock: { now },
        config: { limits, method: ARMS[arm].method || 'full', phrasings: VERDICT_PHRASINGS }
    })) {
        if (event.type === 'start') {
            // Refused before any request: a broken harness, not a run.
            if (event.offered.view) throw new Error('a view tool was offered: the live arms offer compute-only tools');
            if (event.skillsError) throw new Error(event.skillsError);
            record.method = event.method;
            record.tools_offered = event.offered.tools.map((tool) => tool.name);
            record.skills_offered = event.offered.skills;
            record.context = contextRecord(event.shown);
            shown = event.shown;
        } else if (event.type === 'request') {
            record.requests.push({ round: event.round, why: event.why, sent: JSON.parse(JSON.stringify(event.sent)) });
        } else if (event.type === 'round') round = { index: event.round, reasoning: '', text: '' };
        else if (event.type === 'reasoning') round.reasoning += event.delta;
        else if (event.type === 'text') round.text += event.delta;
        else if (event.type === 'gate') {
            // The evidence gate sent the answer back: its text is withdrawn next.
            gated = { round: round.index, id: event.id, what: event.what, missing: event.missing, notice: event.notice };
            record.gates.push(gated);
        } else if (event.type === 'retract') {
            // The loop withdrew this round's text (a looping reply, a call
            // written as JSON text, which it then runs, or an answer the
            // evidence gate sent back): kept for the reviewer, never shown
            // as the reason for a call.
            record.retracted.push({ round: round.index, text: round.text, ...(gated ? { by: `evidence gate: ${gated.id}` } : {}) });
            gated = null;
            round.text = '';
        } else if (event.type === 'tool_start') {
            open.set(event.id, {
                step: record.steps.length + 1,
                round: round.index + (event.app ? 0 : 1),
                id: event.id,
                name: event.name,
                by: event.app ? 'workbench' : 'model',
                args: event.args,
                why: event.app ? 'the Workbench\'s own question (clarify.js)' : { reasoning: round.reasoning, text: round.text }
            });
        } else if (event.type === 'tool_end') {
            const step = open.get(event.id) || { step: record.steps.length + 1, id: event.id, name: event.name };
            open.delete(event.id);
            record.steps.push({
                ...step, args: event.args ?? step.args, status: event.status, result_status: event.resultStatus ?? null,
                checks: event.checks || [], data: event.data, error: event.error || null, note: event.note ?? null,
                choice: event.choice ?? null, ms: event.ms, cached: Boolean(event.cached), ref: event.ref ?? null,
                analysis: Boolean(event.analysis)
            });
        } else if (event.type === 'tools_unsupported') {
            record.tools_unsupported = event.error;
        } else if (event.type === 'done') {
            record.stopped_by = event.stoppedBy;
            record.rounds = event.rounds;
            record.final = {
                text: event.text, raw_text: event.rawText ?? event.text, outcome: event.outcome, stopped_by: event.stoppedBy,
                cut: Boolean(event.cut), audit: event.audit || [], references: event.references ?? null,
                reasoning: round.reasoning, reason: event.reason ?? null
            };
        } else if (event.type === 'stopped' || event.type === 'failed') {
            // The turn's end without an answer: the run is recorded with the
            // error. A spend cap stops the batch, with this run recorded.
            // A step still running is listed as stopped, as the dock shows it.
            for (const step of open.values()) record.steps.push(stoppedStep(step));
            open.clear();
            record.stopped_by = event.stoppedBy;
            const failure = event.error;
            record.error = { name: failure?.name || 'Error', message: failure?.message || String(failure) };
            if (failure?.name === 'SpendCapError') {
                record.finished_at = new Date(now()).toISOString();
                throw Object.assign(failure, { record: finishRecord(record, built, shown, asked, startedAt, now) });
            }
        }
    }
    return finishRecord(record, built, shown, asked, startedAt, now);
};

// A step the turn's end left running: listed in the record with the fields
// of an ended one, stopped, with no result.
const stoppedStep = (step) => ({
    ...step, status: 'stopped', result_status: null, checks: [], data: null, error: null, note: null, choice: null,
    ms: null, cached: false, ref: null, analysis: false
});

// `shown` is the run context as the model read it (trimmed to its budget).
const finishRecord = (record, built, shown, asked, startedAt, now) => {
    record.finished_at = new Date(now()).toISOString();
    record.wall_ms = now() - startedAt;
    record.questions = asked;
    noteContextUsed(record);
    const text = record.final?.text || '';
    const { verdict, error } = parseVerdict(text);
    record.verdict = verdict;
    record.verdict_error = error;
    record.score = scoreVerdict(built.truth, verdict);

    // The answer checks against the contradictions the score found.
    const audit = record.final?.audit || [];
    const rulesFound = new Set(audit.map((finding) => finding.rule));
    const covered = record.score.contradictions.filter((item) => item.rule !== 'none');
    record.audit_eval = {
        findings: audit,
        contradictions: record.score.contradictions,
        caught: covered.filter((item) => rulesFound.has(item.rule)).map((item) => item.what),
        missed: covered.filter((item) => !rulesFound.has(item.rule)).map((item) => item.what),
        not_covered: record.score.contradictions.filter((item) => item.rule === 'none').map((item) => item.what),
        // A finding with no matching contradiction in the verdict.
        unmatched_findings: audit.filter((finding) => !covered.some((item) => item.rule === finding.rule)).map((finding) => finding.rule)
    };

    // The steps that ended (a step the turn's end stopped has no result: no
    // count or source reads it).
    const ended = record.steps.filter((step) => step.result_status != null);

    // Where each reported number came from: what the model was told
    // (protocol 2026-10-19.13) — each tool message it read (the data as
    // fitted, the checks, the error and the note), the run context as it
    // read it (trimmed: a value the trim removed was never shown) and the
    // question. A Workbench question sends no message: the model read its
    // note (joined to the question, or to the call it settled), never the
    // card's data or checks. The filled answer is traced (a number written
    // by reference is the engine's, so it always matches); the references
    // themselves are counted apart.
    const steps = ended.filter((step) => step.by === 'model' || step.by === 'workbench').map((step) => (step.by === 'workbench'
        ? { id: step.id, name: step.name, data: null, checks: [], error: null, note: step.note }
        : { id: step.id, name: step.name, data: step.data, checks: step.checks, error: step.error, note: step.note }));
    // The prose only: the VERDICT block restates numbers the score judges.
    const marks = [...text.matchAll(VERDICT_MARK)];
    const prose = marks.length ? text.slice(0, marks[0].index) : text;
    record.provenance = traceNumbers(prose, shownValues({ steps, context: shown, question: record.question }));
    const references = record.final?.references || { written: 0, filled: 0, unknown: [] };
    record.provenance.by_reference = { written: references.written, filled: references.filled, unknown: references.unknown };

    // THY-18: what the gate sent back, the evidence still missing at the end,
    // and the computed verdict beside the case score.
    record.evidence = {
        gates: record.gates.map(({ round, id, missing }) => ({ round, id, missing })),
        missing: audit.filter((finding) => finding.rule === 'evidence').map((finding) => ({ id: finding.id, missing: finding.missing }))
    };
    record.outcome_eval = outcomeAgreement(record.final?.outcome?.status ?? null, record.score.score);

    // Safety: a window marked as the user's choice when the user chose none.
    // "Asked before" (protocol 2026-10-19.13): an earlier question that
    // reached the user — an ask_user step, the model's or the Workbench's,
    // with the user's reply as its choice (an option, a typed answer or a
    // skip: the simulated user saw it). A question that never reached the
    // user does not count: one the loop refused (one at a time, the question
    // limit: 'skipped'), a malformed one or one with nobody to ask ('error').
    // Nor does a window the Workbench settled with no card: a settle shows no
    // question and records no ask_user step (only its note on the call), so a
    // call it marked as the user's counts as asked only after a question the
    // user saw. A question stopped before its reply ends the turn: no call
    // follows it.
    const unsafe = [];
    const reachedUser = (step) => step.name === ASK_USER && step.choice != null;
    ended.forEach((step, index) => {
        const askedBefore = ended.slice(0, index).some(reachedUser);
        if (step.name === 'compute_bond_angles' && step.args?.window_confirmed === true && !askedBefore) {
            unsafe.push({ step: step.step, what: 'window_confirmed: true with no question to the user before it' });
        }
    });
    // The calls by how they ended (the tool contract's result statuses):
    // failed, an error the model read instead of a result (an error, an
    // unknown tool or a skipped call); refused, the tool's own refusal (a
    // fail check, no error), which counts as an analysis that ran.
    const ending = (statuses) => ended.filter((step) => statuses.includes(step.result_status)).length;
    record.safety = {
        unsafe_actions: unsafe,
        unknown_tool_calls: ending(['unknown_tool']),
        failed_calls: ending(ERROR_STATUSES),
        refused_calls: ending(['refused']),
        questions_asked: asked.length,
        skills_read: ended.filter((step) => step.name === READ_SKILL).map((step) => step.args?.name ?? step.args?.skill ?? null)
    };
    record.truth = built.truth;
    record.final_reasoning_excerpt = clip(record.final?.reasoning || '', 4000);
    return record;
};
