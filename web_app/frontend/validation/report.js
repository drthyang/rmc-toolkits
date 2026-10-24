// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// Every run, side by side: one table per split, a row per case, a column per
// (model, arm), and the summary metrics. Nothing is filtered: a run that
// failed, stopped at the spend cap or wrote no verdict counts as score 0 and
// is listed; the appendix names every run record used. Under each split's
// tables, the comparison beside the case score (comparison.js, plan item 1.6
// part i): the score's parts (the feature-only score among them), every
// component against A+, and every disagreement with the planted truth.
// Before the appendix, what a comparison across models also compares: the
// columns' run policies (runPolicy.js) side by side, and what no policy sets.

import { reproducibility } from './scoring';
import { finalRoundSummary, limitsSummary, SERVER_LABELS, temperatureSummary, windowSummary } from './runPolicy';
import { columnOf, renderComparison } from './comparison';

const mean = (values) => (values.length ? values.reduce((sum, value) => sum + value, 0) / values.length : null);
const fmt = (value, digits = 2) => (value == null ? '–' : value.toFixed(digits));

/** Aggregates run records (live arms and the model-free arms alike). */
export const aggregate = (records) => {
    const groups = new Map();
    records.forEach((record) => {
        const column = columnOf(record);
        const key = `${record.split}|${column}|${record.case}`;
        if (!groups.has(key)) groups.set(key, { split: record.split, column, case: record.case, runs: [] });
        groups.get(key).runs.push(record);
    });
    const cells = [...groups.values()].map((group) => ({
        ...group,
        repro: reproducibility(group.runs.map((run) => ({ verdict: run.verdict, score: run.score }))),
        overclaims: group.runs.filter((run) => run.score?.overclaim).length,
        unsafe: group.runs.reduce((sum, run) => sum + (run.safety?.unsafe_actions?.length || 0), 0),
        caught: group.runs.reduce((sum, run) => sum + (run.audit_eval?.caught?.length || 0), 0),
        covered: group.runs.reduce((sum, run) => sum + (run.audit_eval?.caught?.length || 0) + (run.audit_eval?.missed?.length || 0), 0),
        numbers: group.runs.reduce((sum, run) => sum + (run.provenance?.checkable || 0), 0),
        unsourced: group.runs.reduce((sum, run) => sum + (run.provenance?.unsourced?.length || 0), 0),
        errors: group.runs.filter((run) => run.error).length,
        // THY-18: the gate's send-backs, evidence still missing, numbers by
        // reference, and the computed verdict against the score.
        gates: group.runs.reduce((sum, run) => sum + (run.evidence?.gates?.length || 0), 0),
        missing: group.runs.filter((run) => run.evidence?.missing?.length).length,
        refsFilled: group.runs.reduce((sum, run) => sum + (run.provenance?.by_reference?.filled || 0), 0),
        refsUnknown: group.runs.reduce((sum, run) => sum + (run.provenance?.by_reference?.unknown?.length || 0), 0),
        verdicts: group.runs.map((run) => run.outcome_eval?.agreement ?? 'none')
    }));
    return cells;
};

const distinct = (values) => [...new Set(values)];
const tokens = (count) => Math.round(count).toLocaleString('en-US');

// The judgment subset's situations in the report's order (cases.js
// JUDGMENT_SITUATIONS), with their column headings.
const SITUATION_COLUMNS = Object.freeze([
    ['two_features', 'Two features'],
    ['marginal', 'At the edge of significance'],
    ['two_readings', 'Two tools read differently'],
    ['user_choice', 'The user\'s choice']
]);

/**
 * The judgment subset of one split (plan item 1.6 ii): each column's mean
 * score over the cases of each situation, over the whole subset and over the
 * rest of the split, every run counted (a failed run scores 0), with the
 * number of cases. A case is a judgment case when its record's truth says so
 * (`truth.judgment.situation`, cases.js buildCase). Lines of Markdown; none
 * when the split has no judgment case.
 */
export const renderJudgment = (cells) => {
    const situationOf = (cell) => cell.runs.find((run) => run.truth?.judgment)?.truth.judgment.situation ?? null;
    if (!cells.some(situationOf)) return [];
    const columns = distinct(cells.map((cell) => cell.column)).sort();
    const summary = (mine) => {
        const scores = mine.flatMap((cell) => cell.runs.map((run) => run.score?.score ?? 0));
        return mine.length ? `${fmt(mean(scores))} (${mine.length})` : '–';
    };
    const lines = [
        '**The judgment subset** (plan item 1.6 ii, dev only until its held-out version): cases where A+\'s fixed joining rules are '
            + 'ambiguous or wrong, so a model has room to add. Mean score (cases), every run counted:', '',
        `| Column | ${SITUATION_COLUMNS.map(([, heading]) => heading).join(' | ')} | Judgment subset | The rest |`,
        `| --- | ${SITUATION_COLUMNS.map(() => '---').join(' | ')} | --- | --- |`
    ];
    columns.forEach((column) => {
        const mine = cells.filter((cell) => cell.column === column);
        const cellsOf = (situation) => mine.filter((cell) => situationOf(cell) === situation);
        lines.push(`| ${column} | ${SITUATION_COLUMNS.map(([situation]) => summary(cellsOf(situation))).join(' | ')} | `
            + `${summary(mine.filter(situationOf))} | ${summary(mine.filter((cell) => !situationOf(cell)))} |`);
    });
    lines.push('');
    return lines;
};

// What no run policy sets: each server's or model's own, the same request or not.
export const NOT_IN_THE_POLICY = Object.freeze([
    'sampling beyond the temperature (top_p, top_k, min_p, the repetition penalty): each server\'s or model\'s defaults, never sent',
    'the reasoning (thinking) mode and its effort: each server\'s default',
    'the tokenizer: the guard\'s estimate of a request is the same for every model (characters, not tokens), the real count is not',
    'the weights\' quantization and the server\'s release'
]);

// A column's window over its runs: the num_ctx Ollama was sent (lowest to
// highest), else each run's window (LM Studio's loaded context, before and
// after when it changed; the provider's own).
const columnWindow = (policies) => {
    const kinds = distinct(policies.map((policy) => policy.server.kind ?? 'none'));
    if (kinds.length === 1 && kinds[0] === 'ollama') {
        const used = policies.flatMap((policy) => policy.context.num_ctx_used).filter(Number.isFinite);
        const max = distinct(policies.map((policy) => policy.context.model_max_tokens).filter(Number.isFinite));
        const range = used.length ? (Math.min(...used) === Math.max(...used) ? tokens(used[0]) : `${tokens(Math.min(...used))}–${tokens(Math.max(...used))}`) : 'none sent';
        return `num_ctx ${range}, sized per request${max.length ? ` (the model holds ${max.map(tokens).join(', ')})` : ''}`;
    }
    return distinct(policies.map(windowSummary)).join('; ');
};

/**
 * What a comparison across models also compares (plan item 1.2): each live
 * column's run policy side by side — the server, the window, the reply cap,
 * the final round's tool_choice form, the temperature and the limits — then
 * whether the limits, the temperature and the context guard are the same in
 * every column, and what no policy sets. Lines of Markdown; none without live
 * runs.
 */
export const renderPolicyComparison = (records) => {
    const live = records.filter((record) => record.model);
    if (!live.length) return [];
    const columns = distinct(live.map(columnOf)).sort();
    const lines = ['## What a comparison across models also compares', '',
        'The cases, the question, the prompt, the tools and the scoring are the same in every column. The run policy '
        + '(`record.policy`, validation/runPolicy.js) is not all the same:', '',
        '| Column | Server | Context window | Reply cap | Final round (tool_choice "none") | Temperature | Limits | Policy warnings |',
        '| --- | --- | --- | --- | --- | --- | --- | --- |'];
    columns.forEach((column) => {
        const runs = live.filter((record) => columnOf(record) === column);
        const recorded = runs.filter((record) => record.policy).map((record) => record.policy);
        const missing = runs.length - recorded.length;
        const join = (values) => distinct(values).join('; ') || '–';
        const cells = recorded.length ? [
            join(recorded.map((policy) => (policy.server.kind ? `${SERVER_LABELS[policy.server.kind]}${policy.server.version ? ` ${policy.server.version}` : ''}` : 'not read'))),
            columnWindow(recorded),
            join(recorded.map((policy) => policy.reply_cap?.rule ?? 'server not read')),
            join(recorded.map(finalRoundSummary)),
            join(recorded.map((policy) => temperatureSummary(policy.temperature))),
            join(recorded.map((policy) => limitsSummary(policy.limits))),
            String(recorded.reduce((sum, policy) => sum + (policy.warnings?.length || 0), 0))
        ] : Array(7).fill('not recorded');
        if (missing && recorded.length) cells[0] += ` (${missing} run${missing === 1 ? '' : 's'} with no policy)`;
        lines.push(`| ${column} | ${cells.join(' | ')} |`);
    });
    lines.push('');
    const policies = live.map((record) => record.policy).filter(Boolean);
    const differs = [['the limits', 'limits'], ['the temperature', 'temperature'], ['the context guard', 'context_guard']]
        .filter(([, key]) => distinct(policies.map((policy) => JSON.stringify(policy[key]))).length > 1)
        .map(([name]) => name);
    lines.push(differs.length
        ? `**Not the same in every run: ${differs.join(', ')}.** Read the columns with that in mind.`
        : 'The limits, the temperature and the context guard are the same in every run with a policy.');
    lines.push('');
    lines.push('No run policy sets these, so a comparison across models compares them too:');
    NOT_IN_THE_POLICY.forEach((item) => lines.push(`- ${item};`));
    lines.push('');
    return lines;
};

export const renderReport = (records, { title = 'Copilot validation runs' } = {}) => {
    const cells = aggregate(records);
    const lines = [`# ${title}`, '', `${records.length} run records. Every run is counted; failures score 0.`, ''];
    for (const split of ['dev', 'test']) {
        const inSplit = cells.filter((cell) => cell.split === split);
        if (!inSplit.length) continue;
        const columns = [...new Set(inSplit.map((cell) => cell.column))].sort();
        const cases = [...new Set(inSplit.map((cell) => cell.case))].sort();
        lines.push(`## ${split} set`, '');
        lines.push(`| Case | ${columns.join(' | ')} |`);
        lines.push(`| --- | ${columns.map(() => '---').join(' | ')} |`);
        cases.forEach((id) => {
            const row = columns.map((column) => {
                const cell = inSplit.find((c) => c.case === id && c.column === column);
                if (!cell) return '';
                const flags = [cell.overclaims ? `${cell.overclaims} over-claim` : '', cell.unsafe ? `${cell.unsafe} unsafe` : '', cell.errors ? `${cell.errors} error` : '']
                    .filter(Boolean).join(', ');
                return `${fmt(cell.repro.score_mean)} ± ${fmt(cell.repro.score_sd)} (n=${cell.runs.length}${cell.runs.length > 1 ? `, ${cell.repro.identical ? 'same' : 'varies'}` : ''})${flags ? ` ${flags}` : ''}`;
            });
            lines.push(`| ${id} | ${row.join(' | ')} |`);
        });
        lines.push('');
        lines.push('| Column | Mean score | Null cases over-claimed | Contradictions caught by the checks | Unsafe actions | Reproducible cases | Unsourced numbers |');
        lines.push('| --- | --- | --- | --- | --- | --- | --- |');
        columns.forEach((column) => {
            const mine = inSplit.filter((cell) => cell.column === column);
            const scores = mine.flatMap((cell) => cell.runs.map((run) => run.score?.score ?? 0));
            const nullRuns = mine.flatMap((cell) => cell.runs.filter((run) => run.score?.null_case));
            const caught = mine.reduce((sum, cell) => sum + cell.caught, 0);
            const covered = mine.reduce((sum, cell) => sum + cell.covered, 0);
            const multi = mine.filter((cell) => cell.runs.length > 1);
            lines.push(`| ${column} | ${fmt(mean(scores))} | ${nullRuns.filter((run) => run.score.overclaim).length}/${nullRuns.length} | `
                + `${covered ? `${caught}/${covered}` : '–'} | ${mine.reduce((sum, cell) => sum + cell.unsafe, 0)} | `
                + `${multi.length ? `${multi.filter((cell) => cell.repro.identical).length}/${multi.length}` : '–'} | `
                + `${mine.reduce((sum, cell) => sum + cell.unsourced, 0)}/${mine.reduce((sum, cell) => sum + cell.numbers, 0)} |`);
        });
        lines.push('');
        renderJudgment(inSplit).forEach((line) => lines.push(line));
        lines.push('| Column | Evidence gate: send-backs / runs still missing evidence | Numbers by reference: filled / unknown | Computed verdict: earned / over / under / held |');
        lines.push('| --- | --- | --- | --- |');
        columns.forEach((column) => {
            const mine = inSplit.filter((cell) => cell.column === column);
            const runs = mine.reduce((sum, cell) => sum + cell.runs.length, 0);
            const verdicts = mine.flatMap((cell) => cell.verdicts);
            const count = (kind) => verdicts.filter((value) => value === kind).length;
            lines.push(`| ${column} | ${mine.reduce((sum, cell) => sum + cell.gates, 0)} / ${mine.reduce((sum, cell) => sum + cell.missing, 0)} of ${runs} | `
                + `${mine.reduce((sum, cell) => sum + cell.refsFilled, 0)} / ${mine.reduce((sum, cell) => sum + cell.refsUnknown, 0)} | `
                + `${count('earned')} / ${count('over')} / ${count('under')} / ${count('held')} |`);
        });
        lines.push('');
        renderComparison(records.filter((record) => record.split === split)).forEach((line) => lines.push(line));
    }
    renderPolicyComparison(records).forEach((line) => lines.push(line));
    lines.push('## Every run', '');
    records.map((record) => `- \`${record.id}\`: score ${record.score?.score ?? 0}${record.error ? `, error ${record.error.name}` : ''}`)
        .sort().forEach((line) => lines.push(line));
    lines.push('');
    return lines.join('\n');
};
