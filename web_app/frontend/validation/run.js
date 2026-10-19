// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// The validation's command line (run from web_app/frontend with vite-node, so
// the app's modules load as Vite builds them):
//
//   npx vite-node validation/run.js -- write   --split dev --out <dir>
//   npx vite-node validation/run.js -- arm-a   --split dev [--out <dir>]
//   npx vite-node validation/run.js -- live    --arm B --model qwen3.8-27b --split dev [--runs 3] [--cases id,id] [--out <dir>] [--pca-table]
//   npx vite-node validation/run.js -- report  --in <dir> [--in <dir>…] [--out <file>]
//   npx vite-node validation/run.js -- findings <file.jsonl | folder of .json records>
//
// The held-out test set: `--split test` also needs `--final` and a clean git
// tree (its commit is recorded in every record), so its numbers always come
// from committed code. API models run the test set only, under their spend cap.
//
// `--pca-table` (live, dev set only) gives each run the PCA table, as after a
// visit to the PCA Ellipsoid page, to check that context section's effect;
// its records say so (`view.pcaSites`, `context.sections`) and go to their own
// folder by default (…_pca-table).
//
// Before every live run the run policy's pre-run check (runPolicy.js
// preflight) reads the model's server: a server of another kind, or none, or
// an LM Studio model loaded with less than 48k, stops the batch with what was
// found, what is needed and how to reload the model; nothing is sent and no
// record is written, so the same command resumes once it is fixed.

import { execSync } from 'node:child_process';
import { appendFileSync, existsSync, mkdirSync, readFileSync, readdirSync, statSync, writeFileSync } from 'node:fs';
import { dirname, join, resolve } from 'node:path';
import { fileURLToPath } from 'node:url';
import { providerModel } from '../src/llm/provider/model';
import { buildCase, casesOf } from './cases';
import { caseFiles, loadRun, RUN_NAME } from './runFolder';
import { runArmA, runArmSweep } from './armA';
import { scoreVerdict } from './scoring';
import { truthVerdict } from './verdict';
import { ARMS, PROTOCOL_VERSION, runLiveCase } from './harness';
import { renderTranscript } from './transcript';
import { renderReport } from './report';
import { resolveModel } from './models';
import { createSpendTracker, meteredModel } from './spend';
import { parseFindingFiles, parseFindings } from './findings';
import { noteWindowAfter, preflight, readServer, windowSummary } from './runPolicy';

const HERE = dirname(fileURLToPath(import.meta.url));
const REPO = resolve(HERE, '../../..');

const parseArgs = (argv) => {
    const args = { _: [] };
    for (let i = 0; i < argv.length; i += 1) {
        const token = argv[i];
        if (token === '--') continue;
        if (token.startsWith('--')) {
            const key = token.slice(2);
            const next = argv[i + 1];
            const value = next == null || next.startsWith('--') ? true : (i += 1, next);
            if (key === 'in') args.in = [...(args.in || []), value];
            else args[key] = value;
        } else {
            args._.push(token);
        }
    }
    return args;
};

const commitOf = () => {
    try {
        return execSync('git rev-parse HEAD', { cwd: REPO }).toString().trim();
    } catch {
        return null;
    }
};
const treeClean = () => {
    try {
        return execSync('git status --porcelain', { cwd: REPO }).toString().trim() === '';
    } catch {
        return false;
    }
};

/** The held-out set only from committed code, and only on purpose. */
export const guardSplit = ({ split, final, clean }) => {
    if (split !== 'test') return;
    if (!final) throw new Error('the test set is held out: add --final to run it (reported numbers only, never for tuning)');
    if (!clean) throw new Error('the test set runs from a clean git tree, so every record names the commit it measured');
};

const selectCases = (split, only) => {
    const cases = casesOf(split);
    if (!only) return cases;
    const wanted = String(only).split(',').map((id) => id.trim());
    const unknown = wanted.filter((id) => !cases.some((item) => item.id === id));
    if (unknown.length) throw new Error(`not in the ${split} set: ${unknown.join(', ')}`);
    return cases.filter((item) => wanted.includes(item.id));
};

// The repository's skills, as copilotSkills.js bundles them for the app.
export const loadSkills = () => {
    const root = join(REPO, '.claude/skills');
    const files = {};
    const references = [];
    for (const name of readdirSync(root)) {
        const dir = join(root, name);
        if (!existsSync(join(dir, 'SKILL.md'))) continue;
        for (const file of readdirSync(dir).filter((entry) => entry.endsWith('.md'))) {
            const key = `.claude/skills/${name}/${file}`;
            if (file === 'SKILL.md') files[key] = readFileSync(join(dir, file), 'utf8');
            else references.push(key);
        }
    }
    return { files, references, load: async (path) => readFileSync(join(REPO, path), 'utf8') };
};

/**
 * One live run under the run policy (runPolicy.js): the pre-run check (a
 * RunPolicyError before anything is sent: the server is not the model's
 * kind, or an LM Studio model's loaded window is below the minimum), the run
 * through `base` (the model port: each request as sent carries the context
 * guard's plan, the num_ctx Ollama was sent among it; metered under
 * `tracker` for an API model, as `note`), and the window read again after.
 * A spend cap still throws its SpendCapError, with the run's record.
 */
export const liveRun = async ({
    built, arm, model, runIndex, base, tracker = null, note = null, skills = null, commit = null, pcaTable = false, limits = {}
}) => {
    const server = await preflight(model);
    const port = tracker ? meteredModel(base, tracker, note) : base;
    let record;
    try {
        record = await runLiveCase({ built, arm, model, runIndex, port, skills, commit, pcaTable, limits, server });
    } catch (failure) {
        if (failure?.record) noteWindowAfter(failure.record.policy, await readServer(model));
        throw failure;
    }
    noteWindowAfter(record.policy, await readServer(model));
    return record;
};

const commands = {
    // Each case as <out>/<case>/rmc_run/, the run folder the Workbench opens
    // (Choose folder), named as the harness names every run: the Workbench
    // names a run after its folder, and a case id names the case's kind. The
    // truth sits beside the folder (<out>/<case>/truth.json), not in it.
    write: async (args) => {
        const out = resolve(args.out || 'validation/cases-out');
        for (const item of selectCases(args.split || 'dev', args.cases)) {
            const built = buildCase(item);
            const dir = join(out, built.id);
            const runDir = join(dir, RUN_NAME);
            mkdirSync(runDir, { recursive: true });
            for (const { name, text } of caseFiles(built)) writeFileSync(join(runDir, name), text);
            writeFileSync(join(dir, 'truth.json'), `${JSON.stringify({ id: built.id, split: built.split, question: built.question, truth: built.truth }, null, 1)}\n`);
            console.log(`wrote ${runDir} (truth: ${join(dir, 'truth.json')})`);
        }
    },

    // Arms A and A+ (no model), and the truth (D) as a check of the scoring.
    'arm-a': async (args) => {
        const split = args.split || 'dev';
        guardSplit({ split, final: args.final, clean: split === 'test' ? treeClean() : true });
        const commit = commitOf();
        const records = [];
        for (const item of selectCases(split, args.cases)) {
            const built = buildCase(item);
            const run = await loadRun(built, { ladder: false });
            for (const result of [await runArmA(run), await runArmSweep(run), { arm: 'D', readings: [], verdict: truthVerdict(built.truth) }]) {
                const score = scoreVerdict(built.truth, result.verdict);
                records.push({
                    id: `${built.id}__${result.arm}`, protocol: PROTOCOL_VERSION, commit, case: built.id, split, arm: result.arm,
                    model: null, verdict: result.verdict, score, readings: result.readings, truth: built.truth
                });
                console.log(`${built.id}  arm ${result.arm.padEnd(2)}  score ${score.score.toFixed(2)}  `
                    + `${score.spurious.length ? `spurious: ${score.spurious.map((c) => `${c.kind} ${c.element}`).join(', ')}` : ''}`);
            }
        }
        if (args.out) {
            mkdirSync(dirname(resolve(args.out)), { recursive: true });
            writeFileSync(resolve(args.out), `${JSON.stringify(records, null, 1)}\n`);
        }
    },

    live: async (args) => {
        const arm = args.arm;
        if (!ARMS[arm]) throw new Error('--arm B or --arm C');
        const split = args.split || 'dev';
        const model = resolveModel(args.model, process.env);
        if (model.api) {
            if (!model.splits.includes(split)) throw new Error(`${model.name} runs the ${model.splits.join(', ')} set only`);
            if (!model.arms.includes(arm) && !args['allow-arm']) throw new Error(`${model.name} runs arm ${model.arms.join(', ')} by default (--allow-arm to override within the cap)`);
        }
        guardSplit({ split, final: args.final, clean: split === 'test' ? treeClean() : true });
        const pcaTable = Boolean(args['pca-table']);
        if (pcaTable && split !== 'dev') throw new Error('--pca-table is a dev check: --split dev only');
        const runs = Number(args.runs || model.runsPerCase || 3);
        const out = resolve(args.out || `validation/runs/${new Date().toISOString().slice(0, 10)}_${model.name.replace(/[^\w.-]+/g, '-')}_${arm}_${split}${pcaTable ? '_pca-table' : ''}`);
        mkdirSync(out, { recursive: true });
        const commit = commitOf();
        const skills = arm === 'B' ? loadSkills() : null;
        // The app's production model port, bound to this model (the warmer
        // retry, the context guard, the provider's tool_choice form); one
        // context-size session for the batch.
        const base = providerModel(model, { session: {} });
        let tracker = null;
        if (model.api) {
            const ledgerPath = resolve(args.ledger || 'validation/runs/api-spend.jsonl');
            mkdirSync(dirname(ledgerPath), { recursive: true });
            tracker = createSpendTracker({
                capUSD: model.capUSD,
                prices: model.prices,
                ledger: {
                    load: () => (existsSync(ledgerPath) ? readFileSync(ledgerPath, 'utf8').split('\n').filter(Boolean).map((line) => JSON.parse(line)) : []),
                    append: (entry) => appendFileSync(ledgerPath, `${JSON.stringify(entry)}\n`)
                }
            });
            console.log(`API spend so far $${tracker.spent().toFixed(2)} of $${model.capUSD}`);
        }
        for (const item of selectCases(split, args.cases)) {
            const built = buildCase(item);
            for (let index = 0; index < runs; index += 1) {
                const id = `${built.id}__${arm}__${model.name}__r${index + 1}`.replace(/[^\w.-]+/g, '-');
                if (existsSync(join(out, `${id}.json`))) {
                    console.log(`${id}: already recorded, kept`);
                    continue;
                }
                let record;
                try {
                    record = await liveRun({ built, arm, model, runIndex: index, base, tracker, note: id, skills, commit, pcaTable });
                } catch (failure) {
                    if (failure?.name === 'RunPolicyError') {
                        console.log(`${id}: NOT STARTED (run policy)`);
                        throw failure;
                    }
                    if (failure?.name !== 'SpendCapError') throw failure;
                    record = failure.record;
                    writeFileSync(join(out, `${id}.json`), `${JSON.stringify(record, null, 1)}\n`);
                    writeFileSync(join(out, `${id}.md`), renderTranscript(record));
                    console.log(`${id}: STOPPED, ${failure.message}`);
                    return;
                }
                writeFileSync(join(out, `${id}.json`), `${JSON.stringify(record, null, 1)}\n`);
                writeFileSync(join(out, `${id}.md`), renderTranscript(record));
                console.log(`${id}: score ${record.score.score.toFixed(2)} in ${Math.round(record.wall_ms / 1000)} s, window ${windowSummary(record.policy)}`
                    + `${record.error ? ` (error: ${record.error.message})` : ''}${tracker ? `  [$${tracker.spent().toFixed(2)} spent]` : ''}`);
                record.policy.warnings.forEach((warning) => console.log(`${id}: run policy warning: ${warning}`));
            }
        }
    },

    report: async (args) => {
        const records = [];
        for (const dir of args.in || []) {
            for (const file of readdirSync(resolve(dir)).filter((name) => name.endsWith('.json'))) {
                const content = JSON.parse(readFileSync(join(resolve(dir), file), 'utf8'));
                (Array.isArray(content) ? content : [content]).forEach((record) => records.push(record));
            }
        }
        const text = renderReport(records);
        if (args.out) writeFileSync(resolve(args.out), text);
        else console.log(text);
    },

    // A findings file (JSON lines), or a folder of records, one JSON file each
    // (validation/findings/decisions).
    findings: async (args) => {
        const path = resolve(args._[1]);
        if (statSync(path).isDirectory()) {
            const results = parseFindingFiles(readdirSync(path).map((name) => ({ name, text: readFileSync(join(path, name), 'utf8') })));
            results.forEach((result) => console.log(`${result.file}: ${result.ok ? 'ok' : result.errors.join('; ')}`));
            if (results.some((result) => !result.ok)) process.exitCode = 1;
            return;
        }
        const results = parseFindings(readFileSync(path, 'utf8'));
        results.forEach((result) => console.log(`line ${result.line}: ${result.ok ? 'ok' : result.errors.join('; ')}`));
        if (results.some((result) => !result.ok)) process.exitCode = 1;
    }
};

const main = async () => {
    const args = parseArgs(process.argv.slice(2));
    const command = commands[args._[0]];
    if (!command) {
        console.log('commands: write, arm-a, live, report, findings (see the header of validation/run.js)');
        process.exitCode = 1;
        return;
    }
    await command(args);
};

if (process.env.VITEST !== 'true') {
    main().catch((failure) => {
        console.error(failure?.message || failure);
        process.exitCode = 1;
    });
}
