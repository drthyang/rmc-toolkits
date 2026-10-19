// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// Every run, side by side: one table per split, a row per case, a column per
// (model, arm), and the summary metrics. Nothing is filtered: a run that
// failed, stopped at the spend cap or wrote no verdict counts as score 0 and
// is listed; the appendix names every run record used.

import { reproducibility } from './scoring';

const mean = (values) => (values.length ? values.reduce((sum, value) => sum + value, 0) / values.length : null);
const fmt = (value, digits = 2) => (value == null ? '–' : value.toFixed(digits));

/** Aggregates run records (live arms and the model-free arms alike). */
export const aggregate = (records) => {
    const groups = new Map();
    records.forEach((record) => {
        // A run given the PCA table (run.js live --pca-table, a dev check) is its own column.
        const column = `${record.model?.name ?? 'no model'} · ${record.arm}${record.view?.pcaSites ? ' + PCA table' : ''}`;
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
    }
    lines.push('## Every run', '');
    records.map((record) => `- \`${record.id}\`: score ${record.score?.score ?? 0}${record.error ? `, error ${record.error.name}` : ''}`)
        .sort().forEach((line) => lines.push(line));
    lines.push('');
    return lines.join('\n');
};
