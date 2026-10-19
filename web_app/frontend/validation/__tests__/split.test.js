// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// The held-out test set stays held out: its cases are sealed (only
// holdoutLoader.js may import them, on the --final path or with CI's
// RMC_VALIDATION_HOLDOUT=1; holdoutGate.test.js checks), pinned (a change to
// them, or to the generators they use, fails here until the pin is changed
// in review), they share no id with the dev set, no tool, skill, app source,
// doc or validation text names them or their materials, magnitudes,
// directions and seeds (the holdout scan, holdoutScan.js), the CLI runs them
// only on purpose from committed code, and a supervisor finding on them
// cannot propose tuning.
//
// Without the flag (locally) every test that needs the sealed definitions is
// skipped, and the holdout scan looks for the public ids only; CI runs all.

import { createHash } from 'node:crypto';
import { existsSync, readFileSync, readdirSync, statSync } from 'node:fs';
import { join, relative, resolve } from 'node:path';
import { describe, expect, it } from 'vitest';
import { buildCase, DEV_CASES, JUDGMENT_SITUATIONS, splitOf, TEST_CASE_IDS } from '../cases';
import { parseFindingFiles, parseFindings, validateFinding } from '../findings';
import { ARM_A_TEST_RECORDS, HOLDOUT_DIR, holdoutEnabled, loadTestCases } from '../holdoutLoader';
import { guardSplit, guardTestOutput } from '../run';
import { copilotSkills } from '../../src/copilotSkills';
import { HOLDOUT_ALLOWED, HOLDOUT_LEAKS_PENDING_USER } from './holdoutLists';
import { describeHit, holdoutTokens, readText, scanText } from './holdoutScan';

const FRONTEND = resolve(__dirname, '../..');
const REPO = resolve(FRONTEND, '../..');

// The sealed test cases, only where the loader may open them (CI's flag).
const HOLDOUT = holdoutEnabled();
const TEST_CASES = HOLDOUT ? await loadTestCases() : [];

// sha256 of every test case as built (configuration, partial g(r), truth).
// Change only with a dated note in docs/algorithms/copilot-validation.md:
// results on the old test set are not comparable with the new one.
// 2026-10-19 (protocol 2026-10-19.2): the test breathing case recomposed.
// 2026-10-19 (protocol 2026-10-19.11): five held-out judgment cases join the test set
// (two values then changed after the lead's coincidence check, before any run).
// 2026-10-19 (holdout gate): the definitions sealed, unchanged; the hash holds.
const TEST_SET_SHA256 = '397dd6a2f73e1eea9ff14eca2d8500766743cdeb6daf829e565c1a6ae3d79d50';

// Every file under a folder, never node_modules, dist or the sealed folder
// (which no scan reads).
const SEALED = join(FRONTEND, HOLDOUT_DIR);
const filesUnder = (dir, keep) => readdirSync(dir).flatMap((name) => {
    const path = join(dir, name);
    if (name === 'node_modules' || name === 'dist' || path === SEALED) return [];
    return statSync(path).isDirectory() ? filesUnder(path, keep) : (keep(path) ? [path] : []);
});

describe('the dev / test split', () => {
    it('no shared id: the dev ids and the public test ids, each split\'s own', () => {
        const ids = [...DEV_CASES.map((item) => item.id), ...TEST_CASE_IDS];
        expect(new Set(ids).size).toBe(ids.length);
        expect(DEV_CASES.every((item) => item.id.startsWith('dev-'))).toBe(true);
        expect(TEST_CASE_IDS.every((id) => id.startsWith('test-'))).toBe(true);
        expect(DEV_CASES.every((item) => splitOf(item.id) === 'dev')).toBe(true);
        expect(TEST_CASE_IDS.every((id) => splitOf(id) === 'test')).toBe(true);
        expect(splitOf('test-not-a-case')).toBeNull();
        DEV_CASES.filter((item) => item.devOnly).forEach((item) => expect(typeof item.devOnly, item.id).toBe('string'));
        // Plan item 4.6: the fit-data case.
        expect(DEV_CASES.filter((item) => item.devOnly && !item.judgment).map((item) => item.id)).toEqual(['dev-fit-unconverged']);
    });

    it.skipIf(!HOLDOUT)('the sealed test cases are exactly the public ids, in order', () => {
        expect(TEST_CASES.map((item) => item.id)).toEqual([...TEST_CASE_IDS]);
        TEST_CASES.forEach((item) => expect(buildCase(item, { partials: false }).split, item.id).toBe('test'));
    });

    it.skipIf(!HOLDOUT)('has the same kinds of case on both sides and no shared id', () => {
        const ids = [...DEV_CASES, ...TEST_CASES].map((item) => item.id);
        expect(new Set(ids).size).toBe(ids.length);
        expect(DEV_CASES.every((item) => item.id.startsWith('dev-'))).toBe(true);
        expect(TEST_CASES.every((item) => item.id.startsWith('test-'))).toBe(true);
        // A dev-only case (`devOnly`: why) is for developing tools and
        // skills; it is never a reported number. Outside the judgment subset
        // it has a kind of its own; a judgment case (below) may use the
        // options of a paired kind.
        const paired = DEV_CASES.filter((item) => !item.devOnly);
        expect(paired.map((item) => item.kind).sort()).toEqual(TEST_CASES.map((item) => item.kind).sort());
        expect(paired.length).toBeGreaterThanOrEqual(8);
        expect(TEST_CASES.some((item) => item.devOnly)).toBe(false);
        DEV_CASES.filter((item) => item.devOnly).forEach((item) => {
            expect(typeof item.devOnly, item.id).toBe('string');
            if (!item.judgment) expect(TEST_CASES.some((test) => test.kind === item.kind), item.id).toBe(false);
        });
        // Plan item 4.6: the fit-data case.
        expect(DEV_CASES.filter((item) => item.devOnly && !item.judgment).map((item) => item.id)).toEqual(['dev-fit-unconverged']);
    });

    // Plan item 1.6, part ii: cases where A+'s fixed joining rules are
    // ambiguous or wrong. Since protocol 2026-10-19.11 every dev judgment case
    // has a held-out match, written blind from validation/judgment-spec.md:
    // the same situation, kind and option, one match per dev case.
    it('the judgment subset, its dev side: every situation covered, the truth says which', () => {
        const judged = DEV_CASES.filter((item) => item.judgment);
        expect([...new Set(judged.map((item) => item.situation))].sort()).toEqual(Object.keys(JUDGMENT_SITUATIONS).sort());
        judged.forEach((item) => expect(buildCase(item, { partials: false }).truth.judgment, item.id).toEqual({ situation: item.situation }));
        DEV_CASES.filter((item) => !item.judgment).forEach((item) => {
            expect(item.situation, item.id).toBeUndefined();
            expect(buildCase(item, { partials: false }).truth.judgment, item.id).toBeUndefined();
        });
    });

    it.skipIf(!HOLDOUT)('the judgment subset (plan 1.6 ii): paired by situation, every situation covered, the truth says which', () => {
        const judged = DEV_CASES.filter((item) => item.judgment);
        const matched = TEST_CASES.filter((item) => item.judgment);
        expect(judged.map((item) => item.id)).toEqual([
            'dev-jt-off-centre', 'dev-eight-site', 'dev-octahedral-clusters', 'dev-weak-order-a', 'dev-weak-order-b'
        ]);
        expect(matched.length).toBe(judged.length);
        [...judged, ...matched].forEach((item) => {
            expect(item.judgment, item.id).toBe(true);
            expect(item.devOnly, item.id).toBeUndefined();
            expect(Object.keys(JUDGMENT_SITUATIONS), item.id).toContain(item.situation);
            expect(buildCase(item, { partials: false }).truth.judgment, item.id).toEqual({ situation: item.situation });
        });
        // One match per dev case: the same situation and kind, as many of each.
        const signature = (items) => items.map((item) => `${item.situation}/${item.kind}`).sort();
        expect(signature(matched)).toEqual(signature(judged));
        expect([...new Set(judged.map((item) => item.situation))].sort()).toEqual(Object.keys(JUDGMENT_SITUATIONS).sort());
        // Every other case says nothing of the kind, in its item or its truth.
        [...DEV_CASES, ...TEST_CASES].filter((item) => !item.judgment).forEach((item) => {
            expect(item.situation, item.id).toBeUndefined();
            expect(buildCase(item, { partials: false }).truth.judgment, item.id).toBeUndefined();
        });
        expect(() => buildCase({ ...judged[0], situation: 'nope' }, { partials: false })).toThrow(/unknown judgment situation 'nope'/);
    });

    it.skipIf(!HOLDOUT)('pins the test set: configurations, partials and truths', () => {
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
        // A test-set file written inside the repository goes under the sealed
        // folder or a gitignored output folder.
        expect(() => guardTestOutput('test', join(FRONTEND, ARM_A_TEST_RECORDS))).not.toThrow();
        expect(() => guardTestOutput('test', join(FRONTEND, 'validation/runs/x'))).not.toThrow();
        expect(() => guardTestOutput('test', join(FRONTEND, 'validation/baseline/arm-a-test.json'))).toThrow(/goes under/);
        expect(() => guardTestOutput('test', join(FRONTEND, 'validation/holdout-copy/x.json'))).toThrow(/goes under/);
        expect(() => guardTestOutput('dev', join(FRONTEND, 'validation/baseline/arm-a-dev.json'))).not.toThrow();
    });
});

// The holdout scan. Where it looks: what the tools, the prompts and the app
// say (src: js, jsx, md, and the snapshots of what the model sees), every
// skill file, the Python package, every md under docs/ (the change notes,
// the CHANGELOG, copilot-validation.md and AGENT_PLAN.md among them), the
// top-level docs, and the validation itself outside the sealed folder (its
// js and md files, and its findings: the decision records and the batches).
// The sealed folder is never read. The ids are public, so a file that must
// name a case by its id (copilot-validation.md, cases.js's manifest, a
// decision record on a test case) holds an id-only allowance below; every
// other kind of token is a hit there as anywhere. The two lists of what may
// be named, and where, are holdoutLists.js.
const scannedFiles = () => [
    ...filesUnder(join(FRONTEND, 'src'), (path) => /\.(js|jsx|md)$/.test(path) || path.includes('__snapshots__')),
    ...filesUnder(join(FRONTEND, 'validation'), (path) => /\.(js|jsx|md)$/.test(path) || /[\\/]findings[\\/].*\.jsonl?$/.test(path)),
    ...filesUnder(join(REPO, '.claude/skills'), () => true),
    ...filesUnder(join(REPO, 'rmc_toolkits'), (path) => path.endsWith('.py')),
    ...filesUnder(join(REPO, 'docs'), (path) => path.endsWith('.md')),
    ...readdirSync(REPO).filter((name) => name.endsWith('.md')).map((name) => join(REPO, name))
];

// Every file the model can read: the skill files copilotSkills.js bundles
// for the prompt and read_skill, anything under .claude/skills, and the
// snapshots of what the model is sent.
const SNAPSHOTS = 'web_app/frontend/src/__tests__/modelVisible/__snapshots__/';
const modelVisible = (file) => file.startsWith('.claude/skills/') || file.startsWith(SNAPSHOTS)
    || [...Object.keys(copilotSkills.files), ...copilotSkills.references].includes(file);

const allows = (entry, file, hit) => entry.file === file && entry.kind === hit.kind && (entry.case === hit.case || entry.case === '*');

describe('the holdout scan: no tool, skill, app source, doc or validation text names a test case', () => {
    // Every token with the flag; the public ids alone without it.
    const ID_TOKENS = TEST_CASE_IDS.map((id) => ({ case: id, kind: 'id', name: 'id', text: id }));
    const tokens = HOLDOUT ? TEST_CASES.flatMap((item) => holdoutTokens(item, DEV_CASES)) : ID_TOKENS;
    const active = new Set(tokens.map((token) => token.kind));

    it('scans what it promises: the docs, the change notes, the CHANGELOG, the validation and its records, never the sealed folder', () => {
        const files = scannedFiles().map((path) => relative(REPO, path));
        const must = [
            'docs/algorithms/copilot-validation.md', 'docs/CHANGELOG.md', 'docs/AGENT_PLAN.md', 'AGENTS.md',
            'web_app/frontend/validation/cases.js', 'web_app/frontend/validation/run.js', 'web_app/frontend/validation/judgment-spec.md',
            'web_app/frontend/validation/holdoutLoader.js', 'web_app/frontend/validation/findings/batch-1.jsonl',
            ...readdirSync(join(REPO, 'docs/changes/unreleased')).filter((name) => name.endsWith('.md')).map((name) => `docs/changes/unreleased/${name}`),
            ...readdirSync(join(FRONTEND, 'validation/findings/decisions')).map((name) => `web_app/frontend/validation/findings/decisions/${name}`)
        ];
        must.forEach((file) => expect(files, file).toContain(file));
        const sealed = relative(REPO, SEALED);
        expect(files.filter((file) => file.startsWith(`${sealed}/`))).toEqual([]);
        // An id-only allowance is for ids only.
        HOLDOUT_ALLOWED.filter((entry) => entry.case === '*').forEach((entry) => expect(entry.kind, entry.file).toBe('id'));
    });

    it.skipIf(!HOLDOUT)('derives every test case\'s tokens from its parameters: more than its id', () => {
        expect(tokens.filter((token) => token.kind === 'id')).toEqual(ID_TOKENS);
        TEST_CASES.forEach((item) => {
            const kinds = tokens.filter((token) => token.case === item.id).map((token) => token.kind);
            expect(kinds, item.id).toContain('id');
            expect(kinds.length, item.id).toBeGreaterThanOrEqual(3);
        });
        expect(() => holdoutTokens({ id: 'test-new', kind: 'new_kind', params: {} }, DEV_CASES)).toThrow(/add it to PLANTED/);
    });

    // Plan: the token kinds cover the seeds and every planted magnitude.
    it.skipIf(!HOLDOUT)('every test case yields its seed, and each value and supercell its dev twins do not share', () => {
        TEST_CASES.forEach((item) => {
            const mine = tokens.filter((token) => token.case === item.id);
            expect(mine.filter((token) => token.kind === 'seed').map((token) => token.value), item.id).toEqual([item.params.seed]);
            const twins = DEV_CASES.filter((dev) => !dev.judgment && dev.kind === item.kind);
            const differs = (key) => !twins.some((dev) => JSON.stringify(dev.params[key]) === JSON.stringify(item.params[key]));
            for (const key of ['a', 'delta', 'edge', 'elongation', 'short', 'long', 'bond', 'sigma', 'sigmaCommon', 'sigmaOwn', 'breathe']) {
                if (typeof item.params[key] !== 'number' || !differs(key)) continue;
                expect(mine.some((token) => token.kind === 'magnitude' && token.value === item.params[key]), `${item.id} ${key}`).toBe(true);
            }
            if (item.params.offCentre && differs('offCentre')) {
                expect(mine.some((token) => token.kind === 'magnitude' && token.value === item.params.offCentre.delta), `${item.id} offCentre`).toBe(true);
            }
            if (differs('supercell')) expect(mine.some((token) => token.kind === 'supercell'), `${item.id} supercell`).toBe(true);
        });
    });

    it.skipIf(!HOLDOUT)('the judgment subset removes no test case\'s token (plan 1.6 ii): it is never a kind\'s physics', () => {
        const ordinary = DEV_CASES.filter((item) => !item.judgment);
        expect(ordinary.length).toBeLessThan(DEV_CASES.length);
        expect(TEST_CASES.flatMap((item) => holdoutTokens(item, ordinary))).toEqual(tokens);
    });

    it.skipIf(!HOLDOUT)('finds every test case in its own truth, written as prose', () => {
        TEST_CASES.forEach((item) => {
            const { truth } = buildCase(item, { partials: false });
            const text = [truth.summary, ...truth.traps].join('\n');
            const found = scanText(readText(text), tokens.filter((token) => token.case === item.id && token.kind !== 'id'));
            expect(found.length, item.id).toBeGreaterThan(0);
        });
    });

    // Every scanned file with its hits, each hit with the list entry that covers it.
    let scanned = null;
    const scanAll = () => {
        if (scanned) return scanned;
        scanned = scannedFiles().map((path) => {
            const file = relative(REPO, path);
            const read = readText(readFileSync(path, 'utf8'));
            const hits = scanText(read, tokens).map((hit) => ({
                ...hit,
                text: read.text.slice(...hit.spans[0]),
                leak: HOLDOUT_LEAKS_PENDING_USER.find((entry) => entry.file === file && entry.line === hit.line && entry.case === hit.case && entry.kind === hit.kind),
                allowed: HOLDOUT_ALLOWED.find((entry) => allows(entry, file, hit))
            }));
            return { file, read, hits };
        });
        return scanned;
    };
    const covered = (entry) => scanAll().flatMap(({ hits }) => hits).filter((hit) => hit.leak === entry || (!hit.leak && hit.allowed === entry));

    it('finds nothing outside the two lists, and every entry still matches its count', () => {
        const files = scanAll();
        expect(files.length).toBeGreaterThan(300);
        const unexpected = files.flatMap(({ file, read, hits }) => hits.filter((hit) => !hit.leak && !hit.allowed)
            .map((hit) => describeHit(read, hit, hits, file)));
        expect(unexpected, 'a held-out test case is named (matched text masked as ■): remove it; only text no model can read may be allow-listed').toEqual([]);
        // An entry of a kind the scan looks for here (every kind with the flag, the ids without it).
        [...HOLDOUT_ALLOWED, ...HOLDOUT_LEAKS_PENDING_USER].filter((entry) => active.has(entry.kind)).forEach((entry) => {
            expect(covered(entry).length, `${entry.file}${entry.line ? `:${entry.line}` : ''} ${entry.case} ${entry.kind}: update the count after a review`).toBe(entry.count);
        });
    });

    it('reports every leak pending the user\'s decision, each in text the model can read', () => {
        HOLDOUT_LEAKS_PENDING_USER.forEach((entry) => expect(modelVisible(entry.file), entry.file).toBe(true));
        if (HOLDOUT_LEAKS_PENDING_USER.length) {
            console.warn([
                `HOLDOUT LEAKS PENDING THE USER (${HOLDOUT_LEAKS_PENDING_USER.length}): text the model can read names a held-out test case`,
                ...HOLDOUT_LEAKS_PENDING_USER.map((entry) => `  ${entry.file}:${entry.line}  ${entry.case} (${entry.kind})  ${entry.status}`)
            ].join('\n'));
        }
    });

    it('allows only text no model can read: no skill, bundled file or snapshot, nothing the app can import', () => {
        const appModules = filesUnder(join(FRONTEND, 'src'), (path) => /\.(js|jsx)$/.test(path) && !path.includes('__tests__'))
            .map((path) => readFileSync(path, 'utf8'));
        const snapshots = filesUnder(join(REPO, SNAPSHOTS), () => true).map((path) => readFileSync(path, 'utf8'));
        expect(snapshots.length).toBeGreaterThan(10);
        const bundled = [...Object.keys(copilotSkills.files), ...copilotSkills.references];
        expect(Object.keys(copilotSkills.files).length).toBeGreaterThan(3);
        bundled.forEach((file) => expect(existsSync(join(REPO, file)), `bundled skill file ${file}`).toBe(true));
        HOLDOUT_ALLOWED.forEach((entry) => {
            expect(existsSync(join(REPO, entry.file)), entry.file).toBe(true);
            expect(entry.reason, entry.file).toMatch(/\S.{20,}/);
            expect(modelVisible(entry.file), `${entry.file} is model-visible: a hit there is a leak, not an allowed hit`).toBe(false);
            // No app module imports or reads it (a path ending in its name).
            const stem = entry.file.split('/').pop().replace(/\.[^.]+$/, '');
            const named = new RegExp(`['"\`][^'"\`\\n]*/${stem}(?:\\.[a-z]+)?['"\`]`);
            expect(appModules.some((text) => named.test(text)), `an app module names ${entry.file}`).toBe(false);
            covered(entry).forEach((hit) => {
                const token = new RegExp(`(?<![A-Za-z0-9₀-₉])${hit.text.replace(/[.*+?^${}()|[\]\\]/g, '\\$&')}(?![A-Za-z0-9₀-₉])`);
                expect(snapshots.some((text) => token.test(text)), `${entry.file}: its matched text is in a snapshot of what the model sees`).toBe(false);
            });
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
        // Any held-out id (a test file never writes one: holdoutGate.test.js).
        const id = TEST_CASE_IDS[1];
        const onTest = { ...base, case: id, run: `${id}__B__qwen3.8-27b__r1` };
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

    // Decisions on the validation itself, in the same form, one JSON file per
    // record (validation/findings/decisions/<id>.json), so that parallel PRs
    // never edit the same file. The schema's approval is the user's; a
    // decision taken on the user's delegation says so in its evidence and
    // proposed fix.
    const DECISIONS = join(FRONTEND, 'validation/findings/decisions');
    const decisions = () => {
        const results = parseFindingFiles(readdirSync(DECISIONS).map((name) => ({ name, text: readFileSync(join(DECISIONS, name), 'utf8') })));
        results.forEach((result) => expect(result.errors, result.file).toEqual([]));
        return Object.fromEntries(results.map((result) => [result.id, result.finding]));
    };
    const delegated = (decision) => {
        expect(decision).toMatchObject({ status: 'implemented', needs_user_approval: true, approval: { by: 'user', date: '2026-10-19', batch: 2 } });
        for (const field of ['evidence', 'proposed_fix']) {
            expect(decision[field], field).toMatch(/approved by the manager on the user's delegation, 2026-10-19/i);
        }
    };

    // The one record on a held-out case is named after its case's id, which a
    // test file never writes: it is found by its case.
    const onTestCase = () => Object.entries(decisions()).filter(([, decision]) => splitOf(decision.case) === 'test');

    it('decision records: one JSON file per record, read in id order, each validated by the findings schema', () => {
        const ids = Object.keys(decisions());
        expect([...ids].sort()).toEqual(ids);
        const [[recomposed, { case: heldOut }]] = onTestCase();
        expect(recomposed).toBe(`2026-10-19-${heldOut}-recomposed`);
        expect(ids.filter((id) => id !== recomposed)).toEqual([
            '2026-10-19-aplus-app-rules', '2026-10-19-asked-before', '2026-10-19-engine-verdicts', '2026-10-19-f5-shell-join',
            '2026-10-19-judgment-dev-cases', '2026-10-19-list-sites-split', '2026-10-19-mode-polarisation',
            '2026-10-19-modulation-bar-and-family', '2026-10-19-provenance-told', '2026-10-19-shell-catalogue',
            '2026-10-19-shell-cut-rule'
        ]);
        expect(existsSync(join(FRONTEND, 'validation/findings/decisions.jsonl'))).toBe(false);
        const record = readFileSync(join(DECISIONS, `${recomposed}.json`), 'utf8');
        const read = parseFindingFiles([
            { name: 'README.md', text: 'not a record' },
            { name: '2026-10-20-b.json', text: record },
            { name: '2026-10-20-a.json', text: '{"case":' },
            { name: 'Late Decision.json', text: record }
        ]);
        expect(read.map((result) => result.id)).toEqual(['2026-10-20-a', '2026-10-20-b', 'Late Decision']);
        expect(read.map((result) => result.ok)).toEqual([false, true, false]);
        expect(read[0].errors[0]).toMatch(/^not JSON/);
        expect(read[2].errors).toEqual(['the file name is the record\'s id, YYYY-MM-DD-slug.json (lower-case words and dashes), got Late Decision.json']);
        const extra = parseFindingFiles([{ name: '2026-10-20-c.json', text: JSON.stringify({ ...JSON.parse(record), id: 'x' }) }]);
        expect(extra[0].errors).toEqual(['unknown field \'id\'']);
    });

    it('decisions: the recomposed test breathing case, approved on the user\'s delegation and implemented with its commit', () => {
        const records = onTestCase();
        expect(records).toHaveLength(1);
        const [[, recomposed]] = records;
        expect(TEST_CASE_IDS).toContain(recomposed.case);
        expect(recomposed.kind).toBe('bug');
        delegated(recomposed);
        // Its values are sealed with the cases: the record points there.
        expect(/TEST_SET_SHA256/.test(recomposed.proposed_fix), 'the record points to the sealed cases').toBe(true);
    });

    it('decisions: the engines return their own verdicts (candidate 5), approved on the user\'s delegation, implemented with its tests', () => {
        const verdicts = decisions()['2026-10-19-engine-verdicts'];
        expect(verdicts).toMatchObject({ case: 'dev-jahn-teller', kind: 'bug' });
        delegated(verdicts);
        // Every test it names exists, by its title as written in the source.
        expect(verdicts.tests.length).toBeGreaterThan(0);
        verdicts.tests.forEach((test) => {
            const [path, title] = test.split(' › ');
            expect(readFileSync(join(REPO, path), 'utf8').replace(/\\'/g, '\''), test).toContain(`'${title}'`);
        });
    });

    it('decisions: the judgment dev subset (plan 1.6 ii), approved on the user\'s delegation, implemented with its tests', () => {
        const judgment = decisions()['2026-10-19-judgment-dev-cases'];
        expect(judgment).toMatchObject({ case: 'dev-jt-off-centre', kind: 'bug' });
        expect(DEV_CASES.find((item) => item.id === judgment.case)?.judgment).toBe(true);
        delegated(judgment);
        // Every test it names exists, by its title as written in the source.
        expect(judgment.tests.length).toBeGreaterThan(0);
        judgment.tests.forEach((test) => {
            const [path, title] = test.split(' › ');
            expect(readFileSync(join(REPO, path), 'utf8').replace(/\\'/g, '\''), test).toContain(`'${title}'`);
        });
    });

    it('decisions: list_sites shows the split (plan 4.1), approved on the user\'s delegation, implemented with its known-answer tests', () => {
        const split = decisions()['2026-10-19-list-sites-split'];
        expect(split).toMatchObject({ case: 'dev-offcentre', kind: 'missing-tool' });
        delegated(split);
        // Every test it names exists, by its title as written in the source.
        expect(split.tests.length).toBeGreaterThan(0);
        split.tests.forEach((test) => {
            const [path, title] = test.split(' › ');
            expect(readFileSync(join(REPO, path), 'utf8').replace(/\\'/g, '\''), test).toContain(`'${title}'`);
        });
    });

    it('decisions: arm A+ reads the app\'s rules (candidate 8), approved on the user\'s delegation, implemented with its known-answer tests', () => {
        const rules = decisions()['2026-10-19-aplus-app-rules'];
        expect(rules).toMatchObject({ case: 'dev-jahn-teller', kind: 'bug' });
        delegated(rules);
        expect(rules.tests.length).toBeGreaterThan(0);
        rules.tests.forEach((test) => {
            const [path, title] = test.split(' › ');
            expect(readFileSync(join(REPO, path), 'utf8').replace(/\\'/g, '\''), test).toContain(`'${title}'`);
        });
    });

    it('decisions: the polarisation of an ordered mode (plan 4.3) and the modulation tool\'s bar and family (F2, F3), approved on the user\'s delegation, implemented with their known-answer tests', () => {
        const polarisation = decisions()['2026-10-19-mode-polarisation'];
        expect(polarisation).toMatchObject({ case: 'dev-zone-boundary', kind: 'missing-tool' });
        const barAndFamily = decisions()['2026-10-19-modulation-bar-and-family'];
        expect(barAndFamily).toMatchObject({ case: 'dev-weak-order-b', kind: 'bug' });
        for (const decision of [polarisation, barAndFamily]) {
            delegated(decision);
            expect(decision.tests.length).toBeGreaterThan(0);
            decision.tests.forEach((test) => {
                const [path, title] = test.split(' › ');
                expect(readFileSync(join(REPO, path), 'utf8').replace(/\\'/g, '\''), test).toContain(`'${title}'`);
            });
        }
    });

    it('decisions: the tool transcript\'s stricter readings (candidate 11, D3 and D6), approved on the user\'s delegation, implemented with their known-answer tests', () => {
        const askedBefore = decisions()['2026-10-19-asked-before'];
        const provenanceTold = decisions()['2026-10-19-provenance-told'];
        for (const decision of [askedBefore, provenanceTold]) {
            expect(decision).toMatchObject({ case: 'dev-offcentre', kind: 'bug' });
            delegated(decision);
            expect(decision.tests.length).toBeGreaterThan(0);
            decision.tests.forEach((test) => {
                const [path, title] = test.split(' › ');
                expect(readFileSync(join(REPO, path), 'utf8').replace(/\\'/g, '\''), test).toContain(`'${title}'`);
            });
        }
    });

    it('decisions: one shell reading per pair (candidate 6): the first-peak window approved on the user\'s delegation and implemented; the cut rule only proposed', () => {
        const catalogue = decisions()['2026-10-19-shell-catalogue'];
        expect(catalogue).toMatchObject({ case: 'dev-jahn-teller', kind: 'bug' });
        delegated(catalogue);
        expect(catalogue.tests.length).toBeGreaterThan(0);
        catalogue.tests.forEach((test) => {
            const [path, title] = test.split(' › ');
            expect(readFileSync(join(REPO, path), 'utf8').replace(/\\'/g, '\''), test).toContain(`'${title}'`);
        });
        // The cut-rule recommendation waits for the user: proposed, never approved here.
        const cutRule = decisions()['2026-10-19-shell-cut-rule'];
        expect(cutRule).toMatchObject({ case: 'dev-offcentre', kind: 'orchestration', needs_user_approval: true, status: 'proposed' });
        expect(cutRule.approval).toBeUndefined();
        expect(cutRule.proposed_fix).toMatch(/not approved/i);
    });
});
