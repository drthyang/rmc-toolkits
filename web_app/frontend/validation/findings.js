// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// What a supervisor (a reviewing model or a person reading the transcripts)
// files about a run, in one checked form:
//
//   { case, run, step, kind, evidence, proposed_fix, needs_user_approval }
//
// and, as a finding moves on, where it stands (all optional):
//
//   status: 'proposed' (the default) | 'approved' | 'implemented' | 'rejected'
//   approval: { by: 'user', date: 'YYYY-MM-DD', batch: n }   — who approved it, when, in which batch
//   commit: the commit that implemented it (7–40 hex digits)
//   tests: the regression tests that pin it (['path › test name', …])
//
// - kind: bug (code does the wrong thing), missing-tool (no tool can answer
//   what the question needs), skill-gap (a skill lacks or contradicts a
//   step), orchestration (the agent loop, its limits or guardrails),
//   science-error (a wrong scientific conclusion the tools did not cause).
// - A finding that needs the user's approval is approved or implemented only
//   with its approval; an implemented one names its commit, and an
//   implemented change to code (any kind but skill-gap) names its tests.
// - Bugs are fixed with a regression test and merged on green CI. New tools
//   and skill changes always need the user's approval (batched per cycle).
// - The held-out test set is never used for tuning: a finding on a test case
//   may only report a bug, and it needs the user's approval too.
//
// A batch of findings is a JSON-lines file (findings/batch-1.jsonl,
// parseFindings). Decisions on the validation itself are one JSON file per
// record, named by the record's id (findings/decisions/<id>.json,
// parseFindingFiles), so that parallel changes never edit the same file.

import { splitOf } from './cases';

export const FINDING_KINDS = Object.freeze(['bug', 'missing-tool', 'skill-gap', 'orchestration', 'science-error']);
export const FINDING_STATUSES = Object.freeze(['proposed', 'approved', 'implemented', 'rejected']);
const APPROVAL_REQUIRED = new Set(['missing-tool', 'skill-gap']);
const MIN_EVIDENCE_CHARS = 20;
const ISO_DATE = /^\d{4}-\d{2}-\d{2}$/;
const COMMIT = /^[0-9a-f]{7,40}$/;

// Where a finding stands: status, approval, commit and tests checked together.
const standingErrors = (finding) => {
    const errors = [];
    const status = finding.status ?? 'proposed';
    if (!FINDING_STATUSES.includes(status)) errors.push(`status must be one of ${FINDING_STATUSES.join(', ')}`);
    const { approval } = finding;
    if (approval !== undefined) {
        if (!approval || typeof approval !== 'object' || Array.isArray(approval)) errors.push('approval is { by, date, batch }');
        else {
            if (approval.by !== 'user') errors.push('approval.by is \'user\': only the user approves');
            if (typeof approval.date !== 'string' || !ISO_DATE.test(approval.date)) errors.push('approval.date is YYYY-MM-DD');
            if (!(Number.isInteger(approval.batch) && approval.batch >= 1)) errors.push('approval.batch is the batch number (1, 2, …)');
        }
    }
    if ((status === 'approved' || status === 'implemented') && finding.needs_user_approval === true && approval === undefined) {
        errors.push(`an ${status} finding that needs the user's approval records it (approval)`);
    }
    if (finding.commit !== undefined && (typeof finding.commit !== 'string' || !COMMIT.test(finding.commit))) {
        errors.push('commit is a commit hash (7–40 hex digits)');
    }
    if (status === 'implemented' && finding.commit === undefined) errors.push('an implemented finding names its commit');
    if (finding.tests !== undefined && !(Array.isArray(finding.tests) && finding.tests.every((test) => typeof test === 'string' && test.trim()))) {
        errors.push('tests is a list of test names');
    }
    if (status === 'implemented' && finding.kind !== 'skill-gap' && !(finding.tests?.length)) {
        errors.push('an implemented code change names its regression tests');
    }
    return errors;
};

/** `{ ok, errors }` for one finding. */
export const validateFinding = (finding) => {
    const errors = [];
    if (!finding || typeof finding !== 'object' || Array.isArray(finding)) return { ok: false, errors: ['a finding is an object'] };
    const allowed = ['case', 'run', 'step', 'kind', 'evidence', 'proposed_fix', 'needs_user_approval', 'status', 'approval', 'commit', 'tests'];
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
    errors.push(...standingErrors(finding));
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

// A record's id is its file name without `.json`: the date it was decided
// and a slug, so the id order is the order of the decisions.
export const RECORD_ID = /^\d{4}-\d{2}-\d{2}-[a-z0-9]+(?:-[a-z0-9]+)*$/;

/**
 * Records kept one per file, so that parallel changes never touch the same
 * file (validation/findings/decisions/<id>.json): `files` are `{ name, text }`
 * as read from the folder. Every `*.json` file is one finding, validated by
 * the same schema as a line of a findings file; the others are ignored. In id
 * order; a malformed file or id is an error, never dropped.
 */
export const parseFindingFiles = (files) => (files || [])
    .filter(({ name }) => typeof name === 'string' && name.endsWith('.json'))
    .map(({ name, text }) => ({ id: name.slice(0, -'.json'.length), file: name, text }))
    .sort((a, b) => (a.id < b.id ? -1 : a.id > b.id ? 1 : 0))
    .map(({ id, file, text }) => {
        const idErrors = RECORD_ID.test(id) ? [] : [`the file name is the record's id, YYYY-MM-DD-slug.json (lower-case words and dashes), got ${file}`];
        let finding;
        try {
            finding = JSON.parse(text);
        } catch (failure) {
            return { id, file, finding: null, ok: false, errors: [...idErrors, `not JSON: ${failure.message}`] };
        }
        const { errors } = validateFinding(finding);
        return { id, file, finding, ok: !idErrors.length && !errors.length, errors: [...idErrors, ...errors] };
    });
