// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// The answer corpus (docs/AGENT_PLAN.md item 3.5): real model answers, each
// one file in answerCorpus/entries/, run through the claim reading of the
// evidence gate and the checks on the answer as they are today, on the Demo
// run the answers were given on. Per entry it pins the claims the gate reads,
// the send-back it would make given the steps that ran (or "steps unknown"),
// and the checks that fire — and a person's judgement of each, so a false
// alarm is counted, not guessed. A change to a claim regex that moves a pin
// must update the entry, with a judgement for any new finding.
//
// Fidelity, said by each entry:
//   verbatim     the answer as committed in the repo as the model's own words;
//   transcribed  read off a recording committed in the repo (the README GIF);
//   quoted       a docs or skill quotation of a real answer, possibly shortened
//                by its author: what the whole answer said is not in the repo.
// Every entry is from the Demo run or a dev case of the validation (its live
// runs' records, once they exist); nothing from or about a held-out case.

import { existsSync, readdirSync, readFileSync } from 'node:fs';
import { dirname, join } from 'node:path';
import { fileURLToPath } from 'node:url';
import { beforeAll, describe, expect, it } from 'vitest';
import { DEV_CASES, TEST_CASES } from '../../validation/cases';
import { openCorpusRun, readAnswer, runKeyOf, runSteps, STEPS_UNKNOWN } from './answerCorpus/readAnswer';

const HERE = dirname(fileURLToPath(import.meta.url));
const ENTRIES = join(HERE, 'answerCorpus/entries');
const REPO = join(HERE, '../../../..');
const FIDELITY = ['verbatim', 'transcribed', 'quoted'];

const corpus = readdirSync(ENTRIES).filter((name) => name.endsWith('.json')).sort()
    .map((name) => ({ file: name, entry: JSON.parse(readFileSync(join(ENTRIES, name), 'utf8')) }));

// The findings and send-backs a judgement must cover: "gate:<rule>" for a
// send-back, "audit:<rule>" for each check that fires, "silence" for none.
const judged = (reading) => {
    const on = [
        ...(reading.gateSendBack && reading.gateSendBack !== STEPS_UNKNOWN ? [`gate:${reading.gateSendBack.id}`] : []),
        ...reading.auditFindings.map((rule) => `audit:${rule}`)
    ];
    return on.length ? on : ['silence'];
};

// A JS test file's string literals as their values: concatenations joined,
// escapes read (\\, \', \n), so a verbatim answer committed there is found.
const literalText = (source) => source
    .replace(/'\s*\+\s*'/g, '')
    .replace(/\\([\\'n])/g, (_, char) => (char === 'n' ? '\n' : char));

const valueAt = (object, path) => path.split('.').reduce((value, key) => value?.[key], object);

// Today's false alarms and misses (silence over a conclusion a check or the
// gate is there to read), as `<entry id> <judged on>`. A new one fails here.
const FALSE_ALARMS = [];
const MISSES = ['2026-10-10-gemma4-pooled-no-distortion silence'];

describe('the answer corpus', () => {
    it('holds real answers on the Demo or a dev case, each with its source', () => {
        expect(corpus.length).toBeGreaterThanOrEqual(9);
        const testIds = TEST_CASES.map((item) => item.id);
        const devIds = DEV_CASES.map((item) => item.id);
        for (const { file, entry } of corpus) {
            expect({ file, id: entry.id }).toEqual({ file, id: file.replace(/\.json$/, '') });
            expect(FIDELITY).toContain(entry.fidelity);
            for (const key of ['model', 'date', 'question', 'text', 'stepsSource']) expect(typeof entry[key], `${entry.id}.${key}`).toBe('string');
            expect(/^Demo \(/.test(entry.run) || devIds.includes(entry.run), `${entry.id}: run ${entry.run}`).toBe(true);
            expect(entry.steps === null || Array.isArray(entry.steps), entry.id).toBe(true);
            const raw = JSON.stringify(entry);
            for (const id of testIds) expect(raw.includes(id), `${entry.id} names a held-out case`).toBe(false);
            const source = join(REPO, entry.source.file);
            expect(existsSync(source), `${entry.id}: ${entry.source.file}`).toBe(true);
            if (entry.fidelity === 'quoted') expect(readFileSync(source, 'utf8'), entry.id).toContain(entry.text);
            if (entry.fidelity === 'verbatim') expect(literalText(readFileSync(source, 'utf8')), entry.id).toContain(entry.text);
        }
    });

    describe('read by the current gate and checks', () => {
        const runs = new Map();
        const steps = new Map();
        beforeAll(async () => {
            for (const { entry } of corpus) {
                const key = runKeyOf(entry);
                if (!runs.has(key)) runs.set(key, await openCorpusRun(key));
                if (entry.steps) steps.set(entry.id, await runSteps(runs.get(key).tools, entry.steps));
            }
        }, 120_000);
        const readingOf = (entry) => readAnswer(entry.text, {
            run: runs.get(runKeyOf(entry)),
            ran: steps.get(entry.id)?.ran ?? null,
            facts: steps.get(entry.id)?.facts ?? []
        });

        it.each(corpus.map(({ entry }) => [entry.id, entry]))('%s: the recorded steps run and give the transcript\'s numbers', (id, entry) => {
            if (!entry.steps) return;
            const { results } = steps.get(id);
            results.forEach(({ step, result }) => expect(result?.error ?? null, `${id}: ${step.name}`).toBeNull());
            for (const { step, path, value } of entry.reproduces) {
                expect(valueAt(results[step - 1].result, path), `${id}: step ${step} ${path}`).toEqual(value);
            }
        });

        it.each(corpus.map(({ entry }) => [entry.id, entry]))('%s: claims, send-back and findings as pinned', (id, entry) => {
            expect(readingOf(entry)).toEqual(entry.expected);
        });

        it.each(corpus.map(({ entry }) => [entry.id, entry]))('%s: each finding and send-back is judged', (id, entry) => {
            const on = entry.judgements.map((judgement) => judgement.on);
            expect(on).toEqual(judged(readingOf(entry)));
            for (const judgement of entry.judgements) {
                expect(typeof judgement.right, `${id} ${judgement.on}`).toBe('boolean');
                expect(judgement.reason, `${id} ${judgement.on}`).toMatch(/\S.{10,}/);
            }
        });

        it('opens a dev case too, for answers from the live runs\' records (never a held-out case)', async () => {
            const devRun = await openCorpusRun('dev-thermal-null');
            expect(devRun.offered).toContain('analyze_modulations');
            expect(devRun.context.symmetry?.space_group).toBeTruthy();
            expect(readAnswer('No preferred direction.', { run: devRun, ran: [] })).toEqual({ claims: [], gateSendBack: null, auditFindings: [] });
            await expect(openCorpusRun(TEST_CASES[0].id)).rejects.toThrow(/the Demo or a dev case/);
        }, 120_000);

        it('counts the false alarms and the misses', () => {
            const wrong = (silence) => corpus.flatMap(({ entry }) => entry.judgements
                .filter((judgement) => !judgement.right && (judgement.on === 'silence') === silence)
                .map((judgement) => `${entry.id} ${judgement.on}`));
            expect(wrong(false)).toEqual(FALSE_ALARMS);
            expect(wrong(true)).toEqual(MISSES);
        });
    });
});
