// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// The held-out test set stays held out: its cases are pinned (a change to
// them, or to the generators they use, fails here until the pin is changed
// in review), they share no id with the dev set, no tool, skill, app source
// or doc names them or their materials, magnitudes and directions (the
// holdout scan, holdoutScan.js), the CLI runs them only on purpose from
// committed code, and a supervisor finding on them cannot propose tuning.

import { createHash } from 'node:crypto';
import { readFileSync, readdirSync, statSync } from 'node:fs';
import { join, relative, resolve } from 'node:path';
import { describe, expect, it } from 'vitest';
import { buildCase, DEV_CASES, TEST_CASES } from '../cases';
import { parseFindings, validateFinding } from '../findings';
import { guardSplit } from '../run';
import { describeHit, holdoutTokens, readText, scanText } from './holdoutScan';

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

    it('the CLI runs the test set only with --final, from a clean tree', () => {
        expect(() => guardSplit({ split: 'dev', final: false, clean: false })).not.toThrow();
        expect(() => guardSplit({ split: 'test', final: false, clean: true })).toThrow(/held out: add --final/);
        expect(() => guardSplit({ split: 'test', final: true, clean: false })).toThrow(/clean git tree/);
        expect(() => guardSplit({ split: 'test', final: true, clean: true })).not.toThrow();
    });
});

// The holdout scan. Where it looks: what the tools, the prompts and the app
// say (src: js, jsx, md, and the snapshots of what the model sees), every
// skill file, the Python package, docs/ and the top-level docs. Not scanned:
// the validation itself (web_app/frontend/validation/**, which must name its
// cases) and copilot-validation.md, whose Cases table names every test case.
// docs/AGENT_PLAN.md is scanned: it points at the overlap below without
// naming it.
const NOT_SCANNED = new Set(['docs/algorithms/copilot-validation.md']);
const scannedFiles = () => [
    ...filesUnder(join(FRONTEND, 'src'), (path) => /\.(js|jsx|md)$/.test(path) || path.includes('__snapshots__')),
    ...filesUnder(join(REPO, '.claude/skills'), () => true),
    ...filesUnder(join(REPO, 'rmc_toolkits'), (path) => path.endsWith('.py')),
    ...filesUnder(join(REPO, 'docs'), (path) => path.endsWith('.md')),
    ...readdirSync(REPO).filter((name) => name.endsWith('.md')).map((name) => join(REPO, name))
].filter((path) => !NOT_SCANNED.has(relative(REPO, path)));

// The known hits, each to be reviewed by a person. An entry allows exactly
// `count` hits of one test case's kind of token in one file: a different
// count fails, so a new mention in a listed file is still caught and a
// removed one must leave the list. copilot-validation.md lists these as
// known risks ("The holdout scan"). Nothing here names a test-case value.
const HOLDOUT_ALLOWED = [
    // By design (docs/AGENT_PLAN.md, section A): the test breathing case's
    // cluster is the cluster of the skills' reference compound, whose
    // published distortion this skill describes. A model that has read it may
    // expect a distortion there; the case plants breathing only.
    {
        file: '.claude/skills/coordination-shells/lacunar-spinel-shells.md', case: 'test-breathing', kind: 'cluster', count: 3,
        reason: 'by design: the skills\' reference compound has the test breathing case\'s cluster'
    },
    // FOUND 2026-10-19, to review: the same reference compound's cluster,
    // outside the skills (the 5 K question of the Copilot's real-model check,
    // an engine test's synthetic cluster, the cluster engine's synthetic
    // checks, the Bond Geometry engine's test structures).
    {
        file: 'web_app/frontend/src/llm/README.md', case: 'test-breathing', kind: 'cluster', count: 1,
        reason: 'FOUND 2026-10-19, to review: the reference compound\'s cluster in the real-model check\'s 5 K question'
    },
    {
        file: 'web_app/frontend/src/workers/__tests__/triplets.test.js', case: 'test-breathing', kind: 'cluster', count: 1,
        reason: 'FOUND 2026-10-19, to review: a synthetic cluster of the reference compound in a bond-angle engine test'
    },
    {
        file: 'docs/algorithms/ai-assistant.md', case: 'test-breathing', kind: 'cluster', count: 2,
        reason: 'FOUND 2026-10-19, to review: the cluster engine\'s synthetic checks use the reference compound\'s cluster'
    },
    {
        file: 'docs/algorithms/bond-geometry.md', case: 'test-breathing', kind: 'cluster', count: 1,
        reason: 'FOUND 2026-10-19, to review: the reference compound\'s cluster among the bond-angle engine\'s test geometries'
    },
    // FOUND 2026-10-19, to review: a real-compound fixture of the symmetry
    // finder (an orthorhombic Pnma perovskite) whose composition is the test
    // symmetry-lowering case's; unrelated to the planted shift.
    {
        file: 'web_app/frontend/src/__tests__/fixtures/symmetryStructures.js', case: 'test-symmetry-lowering', kind: 'composition', count: 1,
        reason: 'FOUND 2026-10-19, to review: a symmetry-finder fixture of a real perovskite with the case\'s composition'
    }
];

describe('the holdout scan: no tool, skill, app source or doc names a test case', () => {
    const tokens = TEST_CASES.flatMap((item) => holdoutTokens(item, DEV_CASES));

    it('derives every test case\'s tokens from its parameters: more than its id', () => {
        TEST_CASES.forEach((item) => {
            const kinds = tokens.filter((token) => token.case === item.id).map((token) => token.kind);
            expect(kinds, item.id).toContain('id');
            expect(kinds.length, item.id).toBeGreaterThanOrEqual(3);
        });
        expect(() => holdoutTokens({ id: 'test-new', kind: 'new_kind', params: {} }, DEV_CASES)).toThrow(/add it to PLANTED/);
    });

    it('finds every test case in its own truth, written as prose', () => {
        TEST_CASES.forEach((item) => {
            const { truth } = buildCase(item, { partials: false });
            const text = [truth.summary, ...truth.traps].join('\n');
            const found = scanText(readText(text), tokens.filter((token) => token.case === item.id && token.kind !== 'id'));
            expect(found.length, item.id).toBeGreaterThan(0);
        });
    });

    it('finds nothing outside the reviewed allow-list, and every entry still matches its count', () => {
        const files = scannedFiles();
        expect(files.length).toBeGreaterThan(300);
        const seen = new Map(HOLDOUT_ALLOWED.map((entry) => [entry, 0]));
        const unexpected = [];
        files.forEach((path) => {
            const file = relative(REPO, path);
            const read = readText(readFileSync(path, 'utf8'));
            const hits = scanText(read, tokens);
            hits.forEach((hit) => {
                const entry = HOLDOUT_ALLOWED.find((allowed) => allowed.file === file && allowed.case === hit.case && allowed.kind === hit.kind);
                if (entry) seen.set(entry, seen.get(entry) + 1);
                else unexpected.push(describeHit(read, hit, hits, file));
            });
        });
        expect(unexpected, 'a held-out test case is named (matched text masked as ■); remove it, or add a reviewed HOLDOUT_ALLOWED entry').toEqual([]);
        HOLDOUT_ALLOWED.forEach((entry) => {
            expect(seen.get(entry), `${entry.file} ${entry.case} ${entry.kind}: update the entry's count after a review`).toBe(entry.count);
        });
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
