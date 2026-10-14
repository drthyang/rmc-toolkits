// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// What a supervisor (a reviewing model or a person reading the transcripts)
// files about a run, in one checked form:
//
//   { case, run, step, kind, evidence, proposed_fix, needs_user_approval }
//
// - kind: bug (code does the wrong thing), missing-tool (no tool can answer
//   what the question needs), skill-gap (a skill lacks or contradicts a
//   step), orchestration (the agent loop, its limits or guardrails),
//   science-error (a wrong scientific conclusion the tools did not cause).
// - Bugs are fixed with a regression test and merged on green CI. New tools
//   and skill changes always need the user's approval (batched per cycle).
// - The held-out test set is never used for tuning: a finding on a test case
//   may only report a bug, and it needs the user's approval too.

import { splitOf } from './cases';

export const FINDING_KINDS = Object.freeze(['bug', 'missing-tool', 'skill-gap', 'orchestration', 'science-error']);
const APPROVAL_REQUIRED = new Set(['missing-tool', 'skill-gap']);
const MIN_EVIDENCE_CHARS = 20;

/** `{ ok, errors }` for one finding. */
export const validateFinding = (finding) => {
    const errors = [];
    if (!finding || typeof finding !== 'object' || Array.isArray(finding)) return { ok: false, errors: ['a finding is an object'] };
    const allowed = ['case', 'run', 'step', 'kind', 'evidence', 'proposed_fix', 'needs_user_approval'];
    Object.keys(finding).filter((key) => !allowed.includes(key)).forEach((key) => errors.push(`unknown field '${key}'`));
    const split = typeof finding.case === 'string' ? splitOf(finding.case) : null;
    if (!split) errors.push(`case must be a case id (dev-… or test-…), got ${JSON.stringify(finding.case)}`);
    if (typeof finding.run !== 'string' || !finding.run.trim()) errors.push('run must name the run (its record id)');
    if (!(finding.step === null || (Number.isInteger(finding.step) && finding.step >= 0))) {
        errors.push('step is the transcript step number (0 = the question), or null for the whole run');
    }
    if (!FINDING_KINDS.includes(finding.kind)) errors.push(`kind must be one of ${FINDING_KINDS.join(', ')}`);
    if (typeof finding.evidence !== 'string' || finding.evidence.trim().length < MIN_EVIDENCE_CHARS) {
        errors.push(`evidence must quote the transcript (at least ${MIN_EVIDENCE_CHARS} characters)`);
    }
    if (typeof finding.proposed_fix !== 'string' || !finding.proposed_fix.trim()) errors.push('proposed_fix is required');
    if (typeof finding.needs_user_approval !== 'boolean') errors.push('needs_user_approval must be true or false');
    if (APPROVAL_REQUIRED.has(finding.kind) && finding.needs_user_approval !== true) {
        errors.push(`a ${finding.kind} finding always needs the user's approval`);
    }
    if (split === 'test') {
        if (finding.kind !== 'bug') {
            errors.push(`${finding.case} is a held-out test case: file ${finding.kind} findings from dev cases only (no tuning on the test set)`);
        } else if (finding.needs_user_approval !== true) {
            errors.push('a bug found on a held-out test case needs the user\'s approval before it is fixed');
        }
    }
    return { ok: errors.length === 0, errors };
};

/** Findings from JSON lines, each validated; malformed lines are errors, never dropped. */
export const parseFindings = (text) => String(text || '').split('\n').map((line, index) => ({ line: index + 1, raw: line.trim() }))
    .filter(({ raw }) => raw)
    .map(({ line, raw }) => {
        let finding;
        try {
            finding = JSON.parse(raw);
        } catch (failure) {
            return { line, finding: null, ok: false, errors: [`not JSON: ${failure.message}`] };
        }
        return { line, finding, ...validateFinding(finding) };
    });
