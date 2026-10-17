// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// The held-out test set stays held out: its cases are pinned (a change to
// them, or to the generators they use, fails here until the pin is changed
// in review), they share no id with the dev set, no tool, skill or app
// source names them, the CLI runs them only on purpose from committed code,
// and a supervisor finding on them cannot propose tuning.

import { createHash } from 'node:crypto';
import { readFileSync, readdirSync, statSync } from 'node:fs';
import { join, resolve } from 'node:path';
import { describe, expect, it } from 'vitest';
import { buildCase, DEV_CASES, TEST_CASES } from '../cases';
import { parseFindings, validateFinding } from '../findings';
import { guardSplit } from '../run';

const FRONTEND = resolve(__dirname, '../..');
const REPO = resolve(FRONTEND, '../..');

// sha256 of every test case as built (configuration, partial g(r), truth).
// Change only with a dated note in docs/algorithms/copilot-validation.md:
// results on the old test set are not comparable with the new one.
const TEST_SET_SHA256 = '56dff6c3ad2ae830e7ba1e57bf05fa9bb840eeefc08652f21cd3f09c761b6327';

const filesUnder = (dir, keep) => readdirSync(dir).flatMap((name) => {
    const path = join(dir, name);
    if (name === 'node_modules' || name === 'dist') return [];
    return statSync(path).isDirectory() ? filesUnder(path, keep) : (keep(path) ? [path] : []);
});

describe('the dev / test split', () => {
    it('has the same kinds of case on both sides and no shared id', () => {
        const ids = [...DEV_CASES, ...TEST_CASES].map((item) => item.id);
        expect(new Set(ids).size).toBe(ids.length);
        expect(DEV_CASES.every((item) => item.id.startsWith('dev-'))).toBe(true);
        expect(TEST_CASES.every((item) => item.id.startsWith('test-'))).toBe(true);
        expect(DEV_CASES.map((item) => item.kind).sort()).toEqual(TEST_CASES.map((item) => item.kind).sort());
        expect(DEV_CASES.length).toBeGreaterThanOrEqual(8);
    });

    it('pins the test set: configurations, partials and truths', () => {
        const hash = createHash('sha256');
        TEST_CASES.forEach((item) => {
            const built = buildCase(item);
            hash.update(built.id).update(built.rmc6f).update(built.partials).update(JSON.stringify(built.truth));
        });
        expect(hash.digest('hex')).toBe(TEST_SET_SHA256);
    });

    it('no tool, skill or app source names a test case', () => {
        const ids = TEST_CASES.map((item) => item.id);
        const sources = [
            ...filesUnder(join(FRONTEND, 'src'), (path) => /\.(js|jsx)$/.test(path)),
            ...filesUnder(join(REPO, '.claude/skills'), (path) => path.endsWith('.md')),
            ...filesUnder(join(REPO, 'rmc_toolkits'), (path) => path.endsWith('.py'))
        ];
        expect(sources.length).toBeGreaterThan(50);
        const offenders = sources.filter((path) => ids.some((id) => readFileSync(path, 'utf8').includes(id)));
        expect(offenders).toEqual([]);
    });

    it('the CLI runs the test set only with --final, from a clean tree', () => {
        expect(() => guardSplit({ split: 'dev', final: false, clean: false })).not.toThrow();
        expect(() => guardSplit({ split: 'test', final: false, clean: true })).toThrow(/held out: add --final/);
        expect(() => guardSplit({ split: 'test', final: true, clean: false })).toThrow(/clean git tree/);
        expect(() => guardSplit({ split: 'test', final: true, clean: true })).not.toThrow();
    });
});

describe('supervisor findings', () => {
    const base = {
        case: 'dev-zone-boundary', run: 'dev-zone-boundary__B__qwen3.8-27b__r1', step: 4, kind: 'science-error',
        evidence: 'Step 4: "Ti is randomly off-centred along c" although analyze_modulations reported (0, 0, ½) at > 8σ.',
        proposed_fix: 'Skill: read analyze_modulations before calling a split site disorder.', needs_user_approval: false
    };

    it('accepts a complete dev finding and lists every problem of a bad one', () => {
        expect(validateFinding(base)).toEqual({ ok: true, errors: [] });
        const bad = validateFinding({ case: 'nope', run: '', step: -1, kind: 'typo', evidence: 'short', proposed_fix: '', needs_user_approval: 'yes', extra: 1 });
        expect(bad.ok).toBe(false);
        expect(bad.errors).toHaveLength(8);
    });

    it('new tools and skill changes always need the user\'s approval', () => {
        expect(validateFinding({ ...base, kind: 'skill-gap' }).errors).toEqual(['a skill-gap finding always needs the user\'s approval']);
        expect(validateFinding({ ...base, kind: 'missing-tool', needs_user_approval: true }).ok).toBe(true);
    });

    it('a test case may only report a bug, with approval: no tuning on the test set', () => {
        const onTest = { ...base, case: 'test-zone-boundary', run: 'test-zone-boundary__B__qwen3.8-27b__r1' };
        expect(validateFinding(onTest).errors[0]).toMatch(/held-out test case: file science-error findings from dev cases only/);
        expect(validateFinding({ ...onTest, kind: 'bug' }).errors).toEqual(['a bug found on a held-out test case needs the user\'s approval before it is fixed']);
        expect(validateFinding({ ...onTest, kind: 'bug', needs_user_approval: true }).ok).toBe(true);
    });
});

describe('where a finding stands, and validation batch 1', () => {
    const base = {
        case: 'dev-breathing', run: 'dev-breathing__A+', step: null, kind: 'bug',
        evidence: 'Breathing triangles read "28 % LESS distorted" against independent atoms.',
        proposed_fix: 'A shape null that keeps each cluster\'s breathing.', needs_user_approval: true
    };
    const approval = { by: 'user', date: '2026-10-14', batch: 1 };

    it('an approved or implemented finding records the approval; an implemented one its commit and, for code, its tests', () => {
        expect(validateFinding({ ...base, status: 'approved', approval }).ok).toBe(true);
        expect(validateFinding({ ...base, status: 'approved' }).errors).toEqual(['an approved finding that needs the user\'s approval records it (approval)']);
        expect(validateFinding({ ...base, status: 'implemented', approval }).errors).toEqual([
            'an implemented finding names its commit', 'an implemented code change names its regression tests'
        ]);
        expect(validateFinding({ ...base, status: 'implemented', approval, commit: 'abc1234', tests: ['a.test.js › it works'] }).ok).toBe(true);
        expect(validateFinding({ ...base, kind: 'skill-gap', status: 'implemented', approval, commit: 'abc1234' }).ok).toBe(true);
        expect(validateFinding({ ...base, status: 'done', approval: { by: 'model', date: '14 Oct', batch: 0 }, commit: 'xyz' }).errors).toEqual([
            'status must be one of proposed, approved, implemented, rejected',
            'approval.by is \'user\': only the user approves', 'approval.date is YYYY-MM-DD', 'approval.batch is the batch number (1, 2, …)',
            'commit is a commit hash (7–40 hex digits)'
        ]);
    });

    it('batch 1: every tool change the 2026-10-14 review proposed, approved by the user and implemented, each with its commit', () => {
        const results = parseFindings(readFileSync(join(FRONTEND, 'validation/findings/batch-1.jsonl'), 'utf8'));
        results.forEach((result) => expect(result.errors, `line ${result.line}`).toEqual([]));
        const findings = results.map((result) => result.finding);
        // The tool changes; the skill proposals come with the skills PR (user review).
        const tools = findings.filter((finding) => finding.kind !== 'skill-gap');
        expect(tools).toHaveLength(6);
        findings.forEach((finding) => {
            expect(finding).toMatchObject({ status: 'implemented', needs_user_approval: true, approval: { by: 'user', date: '2026-10-14', batch: 1 } });
            expect(finding.case.startsWith('dev-')).toBe(true);
        });
        // L1–L5 (L4 is two tools).
        expect(tools.map((finding) => finding.kind).sort()).toEqual(['bug', 'bug', 'bug', 'missing-tool', 'missing-tool', 'orchestration']);
        for (const label of ['L1', 'L2', 'L3', 'L4', 'L5']) {
            expect(tools.some((finding) => finding.evidence.includes(label)), label).toBe(true);
        }
    });
});
