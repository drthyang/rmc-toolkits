// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// The computed verdict (THY-18). The Outcome chip under an answer — achieved,
// partial, not verified, not achieved — comes from the deterministic checks,
// never from the model's own grade (which, self-assessed, was "achieved"
// nearly always, even over a cut coordination shell). The model's Outcome line
// only explains; a status word in it is kept as `claimed` and nothing more.
//
// In order, the first that holds:
//   not verified — a check on the answer contradicts it (facts.js), a
//     reference names no value the turn has (references.js), or the answer
//     concludes something whose analysis did not run, or refused (claims.js);
//   not achieved — analyses ran and none of them succeeded;
//   partial — the analyses were cut short (an ending STOPPED_BY marks
//     cutShort, turnPolicy.js: the round budget or the turn's time, or rounds
//     that added nothing), the reply was cut at the length limit, an analysis
//     failed and was not run again successfully, or one warned (use its result
//     with the stated caveat);
//   achieved — otherwise: every analysis the answer rests on ran and passed
//     its checks (or the context answered, and nothing contradicts it).
// The model may lower that, never raise it: when it says its own answer is
// partial or not achieved ("the run cannot tell"), that stands.
// Analyses are the steps that ran (isAnalysis, the one definition: the
// evidence gate's, the Outcome's and the dock footer's, through tool_end's
// `analysis`): not a question to the user, a skill read, a view tool, a call
// to an unknown tool or one skipped by the budget.
// Pure, no React.

import { STOPPED_BY } from '../turnPolicy';

export const VERDICT_STATUSES = ['achieved', 'partial', 'not verified', 'not achieved'];

/**
 * Whether a recorded step is an analysis: it ran (`executed`: not skipped by
 * the budget, not an unknown tool, not a refused question) and it is not a
 * question to the user (`interactive`), a skill read (`reading`) or a view
 * tool (`view`). A failed analysis is one too, and a refused one; the
 * evidence gate also needs it to have answered (no error, not refused).
 */
export const isAnalysis = (step) => Boolean(step?.executed) && !step.interactive && !step.reading && !step.view;

// The findings that make an answer unverified, by rule, as the note says them.
const findingReason = (findings) => {
    const evidence = findings.filter((finding) => finding.rule === 'evidence');
    const references = findings.filter((finding) => finding.rule === 'reference');
    const checks = findings.length - evidence.length - references.length;
    // A refused call backs no claim (claims.js): named as refused, beside what did not run.
    const refused = evidence.some((finding) => finding.refused?.length);
    const named = evidence.flatMap((finding) => (finding.refused?.length ? [...finding.unrun, ...finding.refused] : finding.missing));
    return [
        checks ? `${checks === 1 ? 'a check below contradicts' : `${checks} checks below contradict`} the answer` : null,
        references.length ? 'a number it cites by reference names no value the analyses returned' : null,
        evidence.length ? `the analysis a conclusion needs did not run${refused ? ' or refused' : ''} (${named.join(', ')})` : null
    ].filter(Boolean).join('; ');
};

const listNames = (steps) => [...new Set(steps.map((step) => step.name))].join(', ');

/**
 * The Outcome of one answer. `steps` are the turn's recorded tool steps
 * (`{ name, status: 'ok'|'warn'|'failed', error, executed, interactive, reading, view }`),
 * `findings` the answer's findings (answer checks, unfilled references,
 * missing evidence), `stoppedBy` the turn's ending (a key of STOPPED_BY),
 * `modelOutcome` parseOutcome's reading of the model's line (its status is
 * only kept as `claimed`).
 * Returns `{ status, note, body, claimed, computed: true, reasons }`.
 */
// The model may lower the computed verdict to these, in this order.
const LOWER = ['achieved', 'partial', 'not achieved'];

export const computeOutcome = ({ text = '', steps = [], findings = [], stoppedBy = 'answer', cut = false, modelOutcome = null } = {}) => {
    const reasons = [];
    let status = null;
    if (findings.length) {
        status = 'not verified';
        reasons.push(findingReason(findings));
    }
    const analyses = steps.filter(isAnalysis);
    if (!status && analyses.length && analyses.every((step) => step.status === 'failed' || step.error)) {
        status = 'not achieved';
        reasons.push(`no analysis succeeded (${listNames(analyses)})`);
    }
    if (!status) {
        // A failure counts until the same analysis later succeeds.
        const unresolved = analyses.filter((step, index) => step.status === 'failed'
            && !analyses.slice(index + 1).some((later) => later.name === step.name && later.status !== 'failed'));
        const warned = analyses.filter((step) => step.status === 'warn');
        if (STOPPED_BY[stoppedBy]?.cutShort) reasons.push('the analyses were cut short before the answer');
        if (cut) reasons.push('the reply was cut at the length limit');
        if (unresolved.length) reasons.push(`${listNames(unresolved)} failed`);
        if (warned.length) reasons.push(`${listNames(warned)} warned: its result holds with the stated caveat`);
        if (reasons.length) status = 'partial';
    }
    if (!status) {
        status = 'achieved';
        reasons.push(analyses.length
            ? `${analyses.length} ${analyses.length === 1 ? 'analysis' : 'analyses'} ran and passed ${analyses.length === 1 ? 'its' : 'their'} checks, and no check contradicts the answer`
            : 'answered from the run context, and no check contradicts it');
    }
    // The model's sentence explains; a grade it gave that the checks do not
    // bear out is said, never used — unless it is lower: the model knows when
    // its answer falls short.
    const claimed = modelOutcome?.status ?? null;
    if (claimed && LOWER.includes(status) && LOWER.indexOf(claimed) > LOWER.indexOf(status)) {
        reasons.push(`the model says the answer is ${claimed}`);
        status = claimed;
    }
    const sentence = modelOutcome?.note || '';
    const explanation = claimed && claimed !== status
        ? ` The model said "${claimed}"${sentence ? `: ${sentence}` : '.'}`
        : (sentence ? ` ${sentence}` : '');
    return {
        status,
        note: `${reasons.join('; ')}.${explanation}`,
        body: modelOutcome?.body ?? String(text || '').trim(),
        claimed,
        computed: true,
        reasons
    };
};
