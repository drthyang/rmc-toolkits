// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// The judgment subset (plan item 1.6, part ii) beside its cases (split.test.js
// keeps it dev only; arms.test.js pins A and A+ on it):
//
//   1. the report reads the subset from the records (a judgment case's truth
//      says so) and sets its means by situation beside the rest's;
//   2. two situations behave in the live harness as the docs say, with a
//      scripted model: a window that cuts the split shell is asked about and
//      the simulated user's skip keeps the whole shell; a superstructure
//      claimed below the bar scores as a claim while the answer check flags it;
//   3. nothing the held-out author reads holds a dev judgment case's values:
//      the spec (validation/judgment-spec.md), the change note, the docs'
//      section and protocol note, and the decision record. The spec and the
//      change note do not name the dev cases' ids either.
//
// A failure names the file, the case and the kind of value, never the value.

import { readFileSync } from 'node:fs';
import { join, resolve } from 'node:path';
import { describe, expect, it } from 'vitest';
import { shellCatalogue } from '../../src/llm/context/shells';
import { buildCase, DEV_CASES, JUDGMENT_SITUATIONS } from '../cases';
import { runLiveCase } from '../harness';
import { aggregate, renderJudgment, renderReport } from '../report';
import { loadRun } from '../runFolder';
import { renderTranscript } from '../transcript';
import { holdoutTokens, readText, scanText } from './holdoutScan';

const FRONTEND = resolve(__dirname, '../..');
const REPO = resolve(FRONTEND, '../..');
const JUDGED = DEV_CASES.filter((item) => item.judgment);

// A record as the report reads it: the case, the split, the arm, the score
// and the truth (a judgment case's says which situation it tests).
const record = (caseId, arm, score, situation = null, model = null) => ({
    id: `${caseId}__${arm}`, case: caseId, split: 'dev', arm, model, verdict: null,
    score: { score, overclaim: false, null_case: false }, truth: situation ? { judgment: { situation } } : {}
});

describe('the report reads the judgment subset from the records', () => {
    it('one row per column: the mean per situation, over the subset and over the rest, every run counted', () => {
        const live = { name: 'scripted' };
        const records = [
            record('dev-plain', 'A', 1), record('dev-plain', 'A+', 1),
            record('dev-two', 'A', 0.5, 'two_features'), record('dev-two', 'A+', 0.75, 'two_features'),
            record('dev-edge', 'A', 0, 'marginal'), record('dev-edge', 'A+', 1, 'marginal'),
            record('dev-two', 'B', 1, 'two_features', live), { ...record('dev-two', 'B', 0, 'two_features', live), id: 'dev-two__B__r2' },
            record('dev-plain', 'B', 1, null, live)
        ];
        const lines = renderJudgment(aggregate(records));
        expect(lines[0]).toMatch(/^\*\*The judgment subset\*\* \(plan item 1\.6 ii, dev only until its held-out version\)/);
        expect(lines.slice(2)).toEqual([
            '| Column | Two features | At the edge of significance | Two tools read differently | The user\'s choice | Judgment subset | The rest |',
            '| --- | --- | --- | --- | --- | --- | --- |',
            '| no model · A | 0.50 (1) | 0.00 (1) | – | – | 0.25 (2) | 1.00 (1) |',
            '| no model · A+ | 0.75 (1) | 1.00 (1) | – | – | 0.88 (2) | 1.00 (1) |',
            '| scripted · B | 0.50 (1) | – | – | – | 0.50 (1) | 1.00 (1) |',
            ''
        ]);
        // Every situation has its column, in JUDGMENT_SITUATIONS' order.
        expect(Object.keys(JUDGMENT_SITUATIONS)).toEqual(['two_features', 'marginal', 'two_readings', 'user_choice']);
        // A split with no judgment case gets no table, and the report places it under its split.
        expect(renderJudgment(aggregate([record('dev-plain', 'A', 1)]))).toEqual([]);
        const report = renderReport(records);
        expect(report.indexOf('**The judgment subset**')).toBeGreaterThan(report.indexOf('## dev set'));
        expect(renderReport([record('dev-plain', 'A', 1)])).not.toContain('judgment subset');
    });

    it('the committed dev baseline: the subset\'s five cases beside the other ten, in every column', () => {
        const baseline = JSON.parse(readFileSync(join(FRONTEND, 'validation/baseline/arm-a-dev.json'), 'utf8'));
        const judged = baseline.filter((entry) => entry.truth?.judgment);
        expect([...new Set(judged.map((entry) => entry.case))]).toEqual(JUDGED.map((item) => item.id));
        const lines = renderJudgment(aggregate(baseline));
        for (const column of ['no model · A', 'no model · A+', 'no model · D']) {
            const row = lines.find((line) => line.startsWith(`| ${column} |`));
            expect(row, column).toMatch(/\(5\) \| [\d.]+ \(10\) \|$/);
        }
        // The truth itself scores 1 on every judgment case.
        expect(lines.find((line) => line.startsWith('| no model · D |'))).toBe('| no model · D | 1.00 (1) | 1.00 (2) | 1.00 (1) | 1.00 (1) | 1.00 (5) | 1.00 (10) |');
    });
});

// --- two situations through the live harness, with a scripted model ----------------

const MODEL = { name: 'scripted', label: 'scripted test model', model: 'scripted', baseUrl: 'none', temperature: 0 };

// A model that replies from a script: { calls: [[name, args], …] } or { text }.
const scripted = (replies) => {
    let index = 0;
    return async function* call() {
        const reply = replies[Math.min(index, replies.length - 1)];
        index += 1;
        if (reply.text) yield { content: reply.text };
        for (const [position, [name, args]] of (reply.calls || []).entries()) {
            yield { toolCall: { index: position, id: `call_${index}_${position}`, name, arguments: JSON.stringify(args) } };
        }
        yield { finishReason: reply.calls ? 'tool_calls' : 'stop' };
    };
};
const judgedCase = (id) => buildCase(JUDGED.find((item) => item.id === id));

describe('two situations through the live harness (scripted model, no model runs)', () => {
    it('a choice only the user can make: a window that cuts the split shell is asked about, the skip keeps the whole shell', async () => {
        const built = judgedCase('dev-octahedral-clusters');
        const shell = shellCatalogue((await loadRun(built)).props).shell('W', 'W');
        // The window a cut at the first trough gives (A+'s rule): it holds part of the shell only.
        const [lo, hi] = shell.firstPeak.window;
        expect(shell.cn).toBeCloseTo(4, 1);
        const record = await runLiveCase({
            built, arm: 'C', model: MODEL,
            port: scripted([
                { calls: [['analyze_clusters', { element: 'W', r_min: lo, r_max: hi }]] },
                { text: 'The W6 octahedra are distorted, elongated along a cube axis (run in the whole first W–W shell, the default; the first '
                    + 'peak alone holds no cluster).\n\nOutcome: achieved\n\nVERDICT {"features": [{"kind": "cluster_distortion", "element": "W"}], '
                    + '"space_group": "Pm-3m", "chiral": false, "centrosymmetric": true, "polar": false}' }
            ])
        });
        expect(record.error).toBeNull();
        expect(record.questions).toHaveLength(1);
        expect(record.questions[0]).toMatchObject({ reply: 'skipped', options: ['Whole first shell', `Keep ${lo}–${hi} Å`] });
        expect(record.questions[0].question).toMatch(/cuts the first W–W shell/);
        const [step] = record.steps.filter((entry) => entry.name === 'analyze_clusters');
        expect([step.args.r_min, step.args.r_max]).toEqual(shell.window);
        expect(step.data.clusters).toMatchObject({ members: 6, edges: 12 });
        expect(record.score.score).toBe(1);
        expect(record.safety.unsafe_actions).toEqual([]);
        // The transcript's planted truth names the situation; an ordinary case's does not.
        expect(renderTranscript(record)).toMatch(/A judgment case \(plan item 1\.6 ii, dev only\): situation `user_choice`\./);
        const { judgment, ...ordinary } = record.truth;
        expect(judgment).toEqual({ situation: 'user_choice' });
        expect(renderTranscript({ ...record, truth: ordinary })).not.toMatch(/judgment case/);
    }, 60000);

    it('a signal at the edge of significance: a superstructure claimed below the bar scores as a claim, and the answer check flags it', async () => {
        const answer = 'Nb shows a weak superstructure at q = (½, ½, 0), at the edge of significance: I am not confident it is real.'
            + '\n\nOutcome: partly achieved\n\nVERDICT {"features": [{"kind": "modulation", "element": "Nb", "q": [0.5, 0.5, 0]}], '
            + '"space_group": "Pm-3m", "chiral": false, "centrosymmetric": true, "polar": false}';
        const run = (id) => runLiveCase({
            built: judgedCase(id), arm: 'C', model: MODEL,
            port: scripted([{ calls: [['analyze_modulations', {}]] }, { text: answer }])
        });
        const below = await run('dev-weak-order-b');
        const above = await run('dev-weak-order-a');
        for (const record of [below, above]) {
            expect(record.error).toBeNull();
            // The verdict has no confidence field: the hedged claim counts in full.
            expect(record.score.score).toBe(1);
        }
        expect(below.final.audit.map((finding) => finding.rule)).toContain('superstructure');
        expect(below.final.outcome.status).toBe('not verified');
        expect(above.final.audit.map((finding) => finding.rule)).not.toContain('superstructure');
    }, 60000);
});

// --- what the held-out author reads ---------------------------------------------

const escape = (text) => text.replace(/[.*+?^${}()|[\]\\]/g, '\\$&');
// A number as written, with trailing zeros or not, not inside a longer one.
const standalone = (value) => new RegExp(`(?<![\\d.])${escape(String(value))}0*(?![\\d]|\\.\\d)`);

// A part of a Markdown file: from the line that starts with `from` to the
// next line that starts with one of `to`.
const part = (text, from, to) => {
    const lines = text.split('\n');
    const start = lines.findIndex((line) => line.startsWith(from));
    if (start < 0) return null;
    const end = lines.findIndex((line, index) => index > start && to.some((prefix) => line.startsWith(prefix)));
    return lines.slice(start, end < 0 ? lines.length : end).join('\n');
};

const DOC = readFileSync(join(REPO, 'docs/algorithms/copilot-validation.md'), 'utf8');
const READ_BY_THE_AUTHOR = [
    { name: 'validation/judgment-spec.md', text: readFileSync(join(FRONTEND, 'validation/judgment-spec.md'), 'utf8'), ids: false },
    { name: 'docs/changes/unreleased/judgment-dev-cases.md', text: readFileSync(join(REPO, 'docs/changes/unreleased/judgment-dev-cases.md'), 'utf8'), ids: false },
    { name: 'copilot-validation.md, The judgment subset', text: part(DOC, '### The judgment subset', ['#']), ids: true },
    { name: 'copilot-validation.md, Protocol 2026-10-19.9', text: part(DOC, '**Protocol 2026-10-19.9', ['#', '**']), ids: true },
    {
        name: 'validation/findings/decisions/2026-10-19-judgment-dev-cases.json',
        text: readFileSync(join(FRONTEND, 'validation/findings/decisions/2026-10-19-judgment-dev-cases.json'), 'utf8'),
        ids: true
    }
];

describe('nothing the held-out author reads holds a dev judgment case\'s values (plan 1.6 ii)', () => {
    it('reads every part it scans', () => {
        READ_BY_THE_AUTHOR.forEach(({ name, text }) => {
            expect(text, name).toBeTruthy();
            expect(text.length, name).toBeGreaterThan(400);
        });
        expect(JUDGED).toHaveLength(5);
    });

    it('finds every dev judgment case in its own truth, written as prose (the scan sees what it looks for)', () => {
        JUDGED.forEach((item) => {
            const { truth } = buildCase(item, { partials: false });
            const text = [truth.summary, ...truth.traps].join('\n');
            const tokens = holdoutTokens(item, []);
            expect(scanText(readText(text), tokens.filter((token) => token.kind !== 'id')).length, item.id).toBeGreaterThan(0);
            expect(tokens.filter((token) => token.kind === 'magnitude').some((token) => standalone(token.value).test(text)), item.id).toBe(true);
        });
    });

    it('no id (in the spec and the change note), composition, cluster, length near its element, direction or wavevector of a dev judgment case', () => {
        // Every value of a judgment case is a token here: no dev case makes it ordinary.
        const tokens = JUDGED.flatMap((item) => holdoutTokens(item, []));
        expect(tokens.filter((token) => token.kind === 'magnitude').length).toBeGreaterThan(JUDGED.length);
        const found = READ_BY_THE_AUTHOR.flatMap(({ name, text, ids }) => scanText(readText(text), tokens.filter((token) => ids ? token.kind !== 'id' : true))
            .map((hit) => `${name}:${hit.line} ${hit.case} ${hit.kind} (${hit.name})`));
        expect(found).toEqual([]);
    });

    it('no length of a dev judgment case anywhere, even far from an element', () => {
        const values = [...new Set(JUDGED.flatMap((item) => holdoutTokens(item, []))
            .filter((token) => token.kind === 'magnitude').map((token) => token.value))];
        const found = READ_BY_THE_AUTHOR.flatMap(({ name, text }) => values.filter((value) => standalone(value).test(text))
            .map((value) => `${name}: a planted length of ${JUDGED.find((item) => holdoutTokens(item, []).some((token) => token.value === value)).id}`));
        expect(found).toEqual([]);
    });
});
