// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// The held-out test set stays held out: its cases are pinned (a change to
// them, or to the generators they use, fails here until the pin is changed
// in review), they share no id with the dev set, no tool, skill, app source
// or doc names them or their materials, magnitudes and directions (the
// holdout scan, holdoutScan.js), the CLI runs them only on purpose from
// committed code, and a supervisor finding on them cannot propose tuning.

import { createHash } from 'node:crypto';
import { existsSync, readFileSync, readdirSync, statSync } from 'node:fs';
import { join, relative, resolve } from 'node:path';
import { describe, expect, it } from 'vitest';
import { buildCase, DEV_CASES, JUDGMENT_SITUATIONS, TEST_CASES } from '../cases';
import { parseFindingFiles, parseFindings, validateFinding } from '../findings';
import { guardSplit } from '../run';
import { copilotSkills } from '../../src/copilotSkills';
import { describeHit, holdoutTokens, readText, scanText } from './holdoutScan';

const FRONTEND = resolve(__dirname, '../..');
const REPO = resolve(FRONTEND, '../..');

// sha256 of every test case as built (configuration, partial g(r), truth).
// Change only with a dated note in docs/algorithms/copilot-validation.md:
// results on the old test set are not comparable with the new one.
// 2026-10-19 (protocol 2026-10-19.2): test-breathing recomposed.
// 2026-10-19 (protocol 2026-10-19.11): five held-out judgment cases join the test set
// (two values then changed after the lead's coincidence check, before any run).
const TEST_SET_SHA256 = '397dd6a2f73e1eea9ff14eca2d8500766743cdeb6daf829e565c1a6ae3d79d50';

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
    it('the judgment subset (plan 1.6 ii): paired by situation, every situation covered, the truth says which', () => {
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
// docs/AGENT_PLAN.md is scanned.
const NOT_SCANNED = new Set(['docs/algorithms/copilot-validation.md']);
const scannedFiles = () => [
    ...filesUnder(join(FRONTEND, 'src'), (path) => /\.(js|jsx|md)$/.test(path) || path.includes('__snapshots__')),
    ...filesUnder(join(REPO, '.claude/skills'), () => true),
    ...filesUnder(join(REPO, 'rmc_toolkits'), (path) => path.endsWith('.py')),
    ...filesUnder(join(REPO, 'docs'), (path) => path.endsWith('.md')),
    ...readdirSync(REPO).filter((name) => name.endsWith('.md')).map((name) => join(REPO, name))
].filter((path) => !NOT_SCANNED.has(relative(REPO, path)));

// Two lists, and nothing in either names a test-case value. An entry covers
// exactly `count` hits of one test case's kind of token in one file (and, for
// a leak, on one line): a different count fails, so a new mention there is
// still caught and a removed one must leave the list.
//
// HOLDOUT_ALLOWED: hits the model can never read, in any arm. Model-visible
// text is the system prompt, the tool descriptions and results, read_skill's
// output and every skill file (the SKILL.md files and the reference files
// copilotSkills.js bundles), so an allowed file is none of those, is in no
// snapshot of what the model sees, and no app module can import it (checked
// below). copilot-validation.md lists these as known risks.
const HOLDOUT_ALLOWED = [
    // The four 2026-10-19 entries for the test breathing case's old cluster
    // (the skills' reference compound's) left with protocol 2026-10-19.2: the
    // case was recomposed and nothing names its new cluster.
    //
    // FOUND 2026-10-19, reviewed: a real-compound fixture of the symmetry
    // finder (an orthorhombic Pnma perovskite) whose composition is the test
    // symmetry-lowering case's, unrelated to the planted shift. Imported only by
    // the symmetry tests in src/__tests__, never by the app, the tools or the
    // snapshot tests.
    {
        file: 'web_app/frontend/src/__tests__/fixtures/symmetryStructures.js', case: 'test-symmetry-lowering', kind: 'composition', count: 1,
        reason: 'test-only: a symmetry-finder fixture imported only by src/__tests__/*.test.js; never model-visible'
    }
];

// HOLDOUT_LEAKS_PENDING_USER: hits the model CAN read. They are leaks, not
// allowed hits: printed on every run, pinned line by line so a change shows,
// and left for the user, whose text the skills are. The scan still fails on
// any hit in neither list. copilot-validation.md lists them as open leaks.
// An entry is { file, line, case, kind, count, status }, the status saying
// how the text reaches the model and that it awaits the user's decision.
// Empty since protocol 2026-10-19.2: the three lines of
// coordination-shells/lacunar-spinel-shells.md named the test breathing
// case's old cluster, and the case was recomposed instead of the skill.
const HOLDOUT_LEAKS_PENDING_USER = [];

// Every file the model can read: the skill files copilotSkills.js bundles
// for the prompt and read_skill, anything under .claude/skills, and the
// snapshots of what the model is sent.
const SNAPSHOTS = 'web_app/frontend/src/__tests__/modelVisible/__snapshots__/';
const modelVisible = (file) => file.startsWith('.claude/skills/') || file.startsWith(SNAPSHOTS)
    || [...Object.keys(copilotSkills.files), ...copilotSkills.references].includes(file);

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

    it('the judgment subset removes no test case\'s token (plan 1.6 ii): it is never a kind\'s physics', () => {
        const ordinary = DEV_CASES.filter((item) => !item.judgment);
        expect(ordinary.length).toBeLessThan(DEV_CASES.length);
        expect(TEST_CASES.flatMap((item) => holdoutTokens(item, ordinary))).toEqual(tokens);
    });

    it('finds every test case in its own truth, written as prose', () => {
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
                allowed: HOLDOUT_ALLOWED.find((entry) => entry.file === file && entry.case === hit.case && entry.kind === hit.kind)
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
        [...HOLDOUT_ALLOWED, ...HOLDOUT_LEAKS_PENDING_USER].forEach((entry) => {
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
            expect(modelVisible(entry.file), `${entry.file} is model-visible: a hit there is a leak, not an allowed hit`).toBe(false);
            const stem = entry.file.split('/').pop().replace(/\.[^.]+$/, '');
            const named = new RegExp(`['"\`][^'"\`\\n]*\\b${stem}\\b[^'"\`\\n]*['"\`]`);
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

    it('decision records: one JSON file per record, read in id order, each validated by the findings schema', () => {
        expect(Object.keys(decisions())).toEqual([
            '2026-10-19-aplus-app-rules', '2026-10-19-asked-before', '2026-10-19-engine-verdicts', '2026-10-19-judgment-dev-cases',
            '2026-10-19-list-sites-split', '2026-10-19-mode-polarisation', '2026-10-19-modulation-bar-and-family',
            '2026-10-19-provenance-told', '2026-10-19-shell-catalogue', '2026-10-19-shell-cut-rule', '2026-10-19-test-breathing-recomposed'
        ]);
        expect(existsSync(join(FRONTEND, 'validation/findings/decisions.jsonl'))).toBe(false);
        const record = readFileSync(join(DECISIONS, '2026-10-19-test-breathing-recomposed.json'), 'utf8');
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
        const recomposed = decisions()['2026-10-19-test-breathing-recomposed'];
        expect(recomposed).toMatchObject({ case: 'test-breathing', kind: 'bug' });
        delegated(recomposed);
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
