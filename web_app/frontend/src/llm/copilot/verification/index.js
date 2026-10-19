// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// The verification: one module answers "is this answer verified?" (THY-18;
// architecture candidate #4). One verification serves one turn: the turn
// (turn.js) builds it from the run, the agent loop (agentLoop.js) records
// each tool step the model reads and asks it to judge each answer.
//
//   createVerification({ context, catalog, offered, phrasings, maxSendBacks, method }) → verification
//   verification.record(step, show) → { ref, analysis, ...show(ref) }
//   verification.judge(text, { canSendBack, rounds, stoppedBy, cut }) → { sendBack } | { done } | null
//   verification.claims → the claim rows the turn checks: [{ id, what, prompt, analyses, need }]
//   isAnalysis(step)    → the one definition of an analysis (verdict.js)
//   evidenceSentence(claims), promptGaps(claims), CLAIMS, FACTS   the catalogues (the binding test)
//   shownContext(context), normalizeSymbol(text)                   for the turn and the validation
//
// Configuration:
//   context       the run context (context/runContext.js). {{ctx.…}} resolves
//                 against its shown (trimmed) form only, never a value the
//                 trim removed; the checks read its facts (the watchdog, the
//                 ladder's groups) and the gate its elements;
//   catalog       the run's shell catalogue (context/shells.js); its finite
//                 clusters (clusterElementsOf) need analyze_clusters before
//                 "no order";
//   offered       the tool names the turn offers, tools on or off; null: every
//                 tool; []: no claim row, so no gate, while references are
//                 still filled and the checks and the Outcome computed;
//   phrasings     { [claim id]: { asserted?: RegExp[] (global), literal?: RegExp[] } },
//                 added to the catalogue's own (the validation's verdict
//                 words, validation/verdict.js VERDICT_PHRASINGS);
//   maxSendBacks  send-backs per turn (DEFAULT_MAX_SEND_BACKS);
//   method        the turn's ('off': the send-back names no analysis, claims.js).
//
// record(step, show): once per tool step, in the order the model reads them.
// `step` is the call's ToolResult (the tool contract, tools.js: resultStatus,
// args, error, checks, facts?, …) with the tool's `name` and kind
// (`interactive`, `reading`, `view`); the verification reads the step's
// severity (severityOf: the `status` the Outcome counts) and whether it ran
// (ran: not an unknown tool or a skipped call) from it. `show(ref)` builds the
// tool message — called once, synchronously, returning { data, … } — and the
// ref then resolves to exactly that data (a callback because the ref is
// written into the message before the data is fitted to its budget). Refs are
// r1, r2, … in record order, a step with an error taking none; its facts are
// dropped. A 'refused' step (the tool's own fail check, no error) ran: it is
// an analysis, takes a ref and keeps its facts. `analysis` is
// isAnalysis(step); the gate also needs the step to have no error.
//
// judge(text, …), in a fixed order (the order is the point: a check on the
// unfilled text missed what the filled one says):
//   1. fill the references ({{ctx.…}} from the shown context);
//   2. decide a send-back on the FILLED text: it needs `canSendBack`, an unmet
//      claim (in CLAIMS order) not sent back yet this turn, and fewer than
//      `maxSendBacks` so far; the raw text becomes the standing answer;
//   3. run the answer checks (in FACTS order);
//   4. add the unfilled-reference finding;
//   5. add the missing-evidence findings, tools on, off or unsupported;
//   6. compute the Outcome (verdict.js); the model's Outcome line, read from
//      the filled text, can only lower it.
// `sendBack` is { id, what, missing, nudge, notice }; `done` is the turn's
// `done` event. Blank text judges the standing answer (never sent back
// again), or returns null when there is none.
//
// Errors: createVerification throws a TypeError on a phrasing keyed to an
// unknown claim, or on a non-RegExp or non-global `asserted` pattern (or a
// non-RegExp `literal` one); record throws a TypeError on a step that is no
// tool result (no known resultStatus); record or judge after `done` throws
// (one verification serves one turn). Model text and tool facts never throw.
// Pure, no React.

import { contextToJson } from '../../context/runContext';
import { clusterElementsOf } from '../../context/shells';
import { createReferenceStore, fillReferences } from '../references';
import { TURN_METHODS } from '../method';
import { RESULT_STATUSES, ran as toolRan, severityOf } from '../tools';
import { CLAIM_IDS, claimsOf, evidenceFinding, sendBackOf } from './claims';
import { FACTS, contextAuditFacts } from './facts';
import { parseOutcome } from './outcome';
import { sentencesOf } from './reading';
import { computeOutcome, isAnalysis } from './verdict';

export { CLAIMS, evidenceSentence, promptGaps } from './claims';
export { FACTS } from './facts';
export { normalizeSymbol } from './reading';
export { isAnalysis } from './verdict';

/** Send-backs per turn: an answer goes back at most this often, once per claim. */
export const DEFAULT_MAX_SEND_BACKS = 2;

/**
 * The run context as the model reads it (contextToJson trims it to its
 * budget: a long ladder to its two ends, a long history downsampled): what
 * {{ctx.…}} references name, and the turn's `shown`. Null without a context.
 */
export const shownContext = (context) => (context ? JSON.parse(contextToJson(context)) : null);

const invalid = (message) => {
    throw new TypeError(`createVerification: ${message}`);
};

// The caller's phrasings, checked: known claim ids, RegExp patterns, the
// asserted ones global (they are read with matchAll).
const readPhrasings = (phrasings) => {
    if (phrasings == null) return null;
    if (typeof phrasings !== 'object' || Array.isArray(phrasings)) invalid('phrasings must be { [claim id]: { asserted, literal } } or null');
    Object.entries(phrasings).forEach(([id, entry]) => {
        if (!CLAIM_IDS.includes(id)) invalid(`a phrasing names the unknown claim '${id}' (claims: ${CLAIM_IDS.join(', ')})`);
        if (!entry || typeof entry !== 'object' || Array.isArray(entry)) invalid(`phrasings['${id}'] must be { asserted?, literal? }`);
        const unknown = Object.keys(entry).filter((key) => key !== 'asserted' && key !== 'literal');
        if (unknown.length) invalid(`phrasings['${id}'] has unknown keys: ${unknown.join(', ')}`);
        ['asserted', 'literal'].forEach((kind) => {
            if (entry[kind] == null) return;
            if (!Array.isArray(entry[kind])) invalid(`phrasings['${id}'].${kind} must be an array of RegExp`);
            entry[kind].forEach((pattern) => {
                if (!(pattern instanceof RegExp)) invalid(`phrasings['${id}'].${kind} holds a non-RegExp: ${String(pattern)}`);
                if (kind === 'asserted' && !pattern.global) invalid(`phrasings['${id}'].asserted pattern ${pattern} must be global (g)`);
            });
        });
    });
    return phrasings;
};

/** One turn's verification (see the header). */
export const createVerification = ({
    context = null, catalog = null, offered = null, phrasings = null, maxSendBacks = DEFAULT_MAX_SEND_BACKS, method = 'full'
} = {}) => {
    if (offered != null && !(Array.isArray(offered) && offered.every((name) => typeof name === 'string'))) {
        invalid('offered must be an array of tool names or null');
    }
    if (!(typeof maxSendBacks === 'number' && maxSendBacks >= 0)) invalid('maxSendBacks must be a number, 0 or more');
    if (!TURN_METHODS.includes(method)) invalid(`method must be one of ${TURN_METHODS.join(', ')}`);
    const extra = readPhrasings(phrasings);
    const elements = Object.keys(context?.structure?.elements || {});
    const reader = claimsOf({ offered, extra, elements, values: { clusterElement: clusterElementsOf(catalog) } });
    const references = createReferenceStore({ context: shownContext(context) });
    const auditContext = contextAuditFacts(context);
    // The turn so far: its recorded steps, its tools' facts, the claims sent
    // back, the gate's events for `done`, the answer sent back last.
    const steps = [];
    const facts = [];
    const sentBack = new Set();
    const gates = [];
    let standing = null;
    let judged = false;

    const open = (operation) => {
        if (judged) throw new Error(`verification.${operation}: this turn's answer is judged; one verification serves one turn`);
    };
    // The calls the gate counts: analyses that ran without an error.
    const ran = () => steps.filter((step) => isAnalysis(step) && !step.error).map((step) => ({ name: step.name, args: step.args }));
    const checks = (answer) => {
        if (!answer.trim()) return [];
        const reading = { text: answer, sentences: sentencesOf(answer) };
        return FACTS.flatMap((fact) => (fact.contradicts ? fact.contradicts(reading, facts, auditContext) : []));
    };
    const finish = (text, { rounds, stoppedBy, cut }) => {
        const filled = fillReferences(text, references);
        const answer = filled.text;
        const findings = [
            ...checks(answer),
            ...(filled.unknown.length ? [{
                rule: 'reference',
                text: `The answer cites ${filled.unknown.map((name) => `{{${name}}}`).join(', ')}, which `
                    + `${filled.unknown.length === 1 ? 'names' : 'name'} no value this question's analyses returned: those numbers are missing.`
            }] : []),
            // With tools off (or unsupported) too: a conclusion no analysis
            // backs is not verified, whether or not one could have run.
            ...reader.unmet(answer, ran()).map(evidenceFinding)
        ];
        const outcome = computeOutcome({ text: answer, steps, findings, stoppedBy, cut, modelOutcome: parseOutcome(answer) });
        judged = true;
        return {
            type: 'done',
            text: answer,
            ...(answer !== text ? { rawText: text } : {}),
            outcome,
            rounds,
            stoppedBy,
            cut,
            audit: findings,
            references: { written: filled.written, filled: filled.filled, unknown: filled.unknown },
            gates: [...gates]
        };
    };

    return Object.freeze({
        claims: Object.freeze(reader.rows.map((claim) => Object.freeze({
            id: claim.id,
            what: claim.what,
            prompt: claim.prompt,
            analyses: Object.freeze([...new Set(claim.backedBy.requirements.map((item) => item.tool))]
                .filter((tool) => !offered || offered.includes(tool))),
            need: claim.backedBy.need
        }))),

        record(step, show) {
            open('record');
            if (!step || typeof step !== 'object' || typeof step.name !== 'string') throw new TypeError('verification.record: step must be { name, … }');
            if (!RESULT_STATUSES.includes(step.resultStatus)) {
                throw new TypeError(`verification.record: step must be a tool result (resultStatus one of ${RESULT_STATUSES.join(', ')}; tools.js toolResult), got ${JSON.stringify(step.resultStatus)}`);
            }
            if (typeof show !== 'function') throw new TypeError('verification.record: show must be a function of the ref');
            const ref = step.error ? null : `r${references.size() + 1}`;
            const shown = show(ref);
            if (!shown || typeof shown !== 'object') throw new TypeError('verification.record: show(ref) must return { data, … }');
            if (ref) references.put(ref, shown.data);
            // A failed step's facts are dropped; facts that are no object are no facts.
            if (!step.error && step.facts && typeof step.facts === 'object') facts.push(step.facts);
            const entry = {
                name: step.name, args: step.args, status: severityOf(step), error: step.error || null,
                interactive: Boolean(step.interactive), reading: Boolean(step.reading), view: Boolean(step.view),
                executed: toolRan(step), ref
            };
            steps.push(entry);
            return { ref, analysis: isAnalysis(entry), ...shown };
        },

        judge(text, { canSendBack = false, rounds = 0, stoppedBy = 'answer', cut = false } = {}) {
            open('judge');
            const raw = typeof text === 'string' ? text : String(text ?? '');
            const options = { rounds, stoppedBy, cut };
            // No answer: the one the gate sent back stands, if any.
            if (!raw.trim()) return standing == null ? null : { done: finish(standing, options) };
            if (canSendBack && gates.length < maxSendBacks) {
                // On the answer as it will read: references filled.
                const claim = reader.unmet(fillReferences(raw, references).text, ran()).find((item) => !sentBack.has(item.id));
                if (claim) {
                    sentBack.add(claim.id);
                    gates.push({ round: rounds, id: claim.id, missing: claim.missing });
                    standing = raw;
                    return { sendBack: sendBackOf(claim, method) };
                }
            }
            return { done: finish(raw, options) };
        }
    });
};
