// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// A run record as a Markdown transcript a reviewer can read offline, with
// nothing else open: the case, what the model was asked, every step (who
// decided it, why, the arguments, the checks and the result it read), the
// answer, its verdict and score, the answer checks, where every number came
// from, safety, and last, the planted truth the model never saw. Steps are
// numbered: a finding's `step` is that number (0 = the question). The header
// also states the run policy the run ran under (runPolicy.js).

import { policyLines } from './runPolicy';

const fence = (value) => `\`\`\`json\n${JSON.stringify(value, null, 1)}\n\`\`\``;
const quote = (text) => String(text || '').trim().split('\n').map((line) => `> ${line}`).join('\n');
const clip = (text, limit) => (text && text.length > limit ? `${text.slice(0, limit)}\n… [${text.length - limit} more characters in the JSON record]` : text);

const scoreTable = (score) => {
    const rows = ['| Item | Result |', '| --- | --- |'];
    (score.features || []).forEach((feature) => {
        const parts = [feature.found ? 'found' : 'missed'];
        if (feature.direction_error_deg != null) parts.push(`direction off by ${feature.direction_error_deg.toFixed(1)}°${feature.direction_ok ? '' : ' (over 15°)'}`);
        else if (feature.found && 'direction_ok' in feature) parts.push('no direction given');
        if ('magnitude_ok' in feature) parts.push(`magnitude ${feature.magnitude_claimed_A ?? 'not given'}${feature.magnitude_ok ? ' ok' : ' off'}`);
        if ('q_ok' in feature) parts.push(`q ${JSON.stringify(feature.q_claimed)}${feature.q_ok ? ' ok' : ' wrong'}`);
        rows.push(`| planted ${feature.planted || feature.kind} (${feature.element}) | ${parts.join('; ')} |`);
    });
    rows.push(`| spurious claims | ${score.spurious?.length ? score.spurious.map((c) => `${c.kind} (${c.element})`).join(', ') : 'none'} |`);
    if (score.tolerated?.length) rows.push(`| tolerated (true, partial) | ${score.tolerated.map((c) => `${c.kind} (${c.element})`).join(', ')} |`);
    Object.entries(score.symmetry || {}).forEach(([field, result]) => rows.push(`| ${field} | ${result} |`));
    rows.push(`| **case score** | **${score.score}** ${score.format_ok ? '' : '(no parseable VERDICT)'} |`);
    return rows.join('\n');
};

export const renderTranscript = (record) => {
    const lines = [];
    lines.push(`# ${record.id}`);
    lines.push('');
    lines.push(`- Case: \`${record.case}\` (${record.split} set) · Arm ${record.arm}: ${record.arm_label}`);
    lines.push(`- Model: ${record.model.label} (\`${record.model.id}\`, temperature ${record.model.temperature}) · run ${record.run_index}`);
    lines.push(`- Protocol ${record.protocol} · commit ${record.commit || 'unknown'} · ${record.started_at} → ${record.finished_at} (${Math.round((record.wall_ms || 0) / 1000)} s)`);
    // The rules the run ran under (runPolicy.js).
    policyLines(record.policy).forEach((line) => lines.push(line));
    lines.push(`- Tools offered: ${record.tools_offered.join(', ')}`);
    lines.push(`- Skills offered: ${record.skills_offered.length ? record.skills_offered.join(', ') : 'none (ablation)'}`);
    if (record.view && record.context) {
        const view = record.view;
        lines.push(`- View: Detected SG at ${view.toleranceA} Å${view.ladder ? ' with its ladder' : ''}, `
            + `${view.hiddenPlots.length ? `${view.hiddenPlots.length} plot(s) hidden` : 'nothing hidden'}, `
            + `${view.pcaSites ? 'the PCA table' : 'no PCA table'}, Live Data ${view.liveData ? 'on' : 'off'}`);
        lines.push(`- Context sections: ${record.context.sections.join(', ') || 'none'}`
            + `${record.context.absent.length ? ` (absent: ${record.context.absent.join(', ')})` : ''}`
            + ` · ${record.context.earlier_turns ? `${record.context.earlier_turns} earlier turns` : 'one fresh question, no earlier turns'}`);
    }
    if (record.error) lines.push(`- **Run error:** ${record.error.name}: ${record.error.message}`);
    lines.push('');
    lines.push('## Step 0 · the question');
    lines.push('');
    lines.push(quote(record.question));
    record.steps.forEach((step) => {
        lines.push('');
        lines.push(`## Step ${step.step} · ${step.name}${step.by === 'workbench' ? ' (asked by the Workbench)' : ''} · round ${step.round ?? '?'}`);
        lines.push('');
        if (step.by === 'model' && step.why) {
            if (step.why.reasoning?.trim()) lines.push(`**Model reasoning before the call:**\n\n${quote(clip(step.why.reasoning, 2500))}\n`);
            if (step.why.text?.trim()) lines.push(`**Model text before the call:**\n\n${quote(clip(step.why.text, 1500))}\n`);
            if (!step.why.reasoning?.trim() && !step.why.text?.trim()) lines.push('*No reasoning or text was streamed before this call.*\n');
        } else if (typeof step.why === 'string') {
            lines.push(`*${step.why}*\n`);
        }
        lines.push(`**Arguments:** \`${JSON.stringify(step.args)}\` · status ${step.status}${step.cached ? ' (cached repeat)' : ''} · ${step.ms ?? 0} ms`);
        if (step.error) lines.push(`\n**Error:** ${step.error}`);
        if (step.checks?.length) {
            lines.push('\n**Checks the model read:**\n');
            step.checks.forEach((check) => lines.push(`- ${check.level}: ${check.text}`));
        }
        if (step.data != null) lines.push(`\n<details><summary>Result the model read</summary>\n\n${fence(step.data)}\n\n</details>`);
    });
    lines.push('');
    lines.push('## The answer');
    lines.push('');
    (record.retracted || []).forEach((item) => {
        lines.push(item.by
            ? `*Round ${item.round}: the ${item.by} sent the answer below back to the model:*\n`
            : `*Round ${item.round}: the loop withdrew the text below (a looping reply, or a call written as text):*\n`);
        lines.push(`${quote(clip(item.text || '(empty)', 1500))}\n`);
    });
    (record.gates || []).forEach((gate) => lines.push(`- Evidence gate, round ${gate.round}: ${gate.notice}`));
    if (record.gates?.length) lines.push('');
    if (record.final) {
        const outcome = record.final.outcome;
        lines.push(`Stopped by: ${record.final.stopped_by}${record.final.reason ? ` (${record.final.reason})` : ''}`
            + `${record.final.cut ? ' (cut at the length limit)' : ''} · rounds: ${record.rounds} · outcome (computed): ${outcome?.status ?? 'none'}`
            + `${outcome?.claimed && outcome.claimed !== outcome.status ? ` (the model said "${outcome.claimed}")` : ''}`
            + `${record.outcome_eval ? ` · against the score: ${record.outcome_eval.agreement}` : ''}`);
        if (outcome?.note) lines.push(`\nWhy: ${outcome.note}`);
        lines.push('');
        lines.push(quote(record.final.text || '(empty)'));
        const refs = record.final.references;
        if (refs?.written) {
            lines.push('');
            lines.push(`Numbers by reference: ${refs.filled} of ${refs.written} filled from the analyses`
                + `${refs.unknown.length ? `; not filled: ${refs.unknown.map((name) => `{{${name}}}`).join(', ')}` : ''}.`);
            if (record.final.raw_text && record.final.raw_text !== record.final.text) {
                lines.push(`\n<details><summary>The answer as the model wrote it</summary>\n\n${quote(clip(record.final.raw_text, 4000))}\n\n</details>`);
            }
        }
    } else {
        lines.push('*No answer.*');
    }
    lines.push('');
    lines.push('## Verdict and score');
    lines.push('');
    lines.push(record.verdict ? fence(record.verdict) : `*No verdict: ${record.verdict_error}.*`);
    lines.push('');
    lines.push(scoreTable(record.score));
    lines.push('');
    lines.push('## Answer checks');
    lines.push('');
    const audit = record.audit_eval;
    lines.push(`- Findings under the answer: ${audit.findings.length ? audit.findings.map((f) => `${f.rule}: ${f.text}`).join(' · ') : 'none'}`);
    lines.push(`- Contradictions in the verdict a check covers: caught ${audit.caught.length ? audit.caught.join(', ') : 'none'}; missed ${audit.missed.length ? audit.missed.join(', ') : 'none'}`);
    if (audit.not_covered.length) lines.push(`- Contradictions no check covers: ${audit.not_covered.join(', ')}`);
    if (audit.unmatched_findings.length) lines.push(`- Findings with no matching contradiction in the verdict: ${audit.unmatched_findings.join(', ')}`);
    lines.push('');
    lines.push('## Numbers in the answer');
    lines.push('');
    const prov = record.provenance;
    lines.push(`${prov.sourced} of ${prov.checkable} checkable numbers match a value the model was shown.`);
    if (prov.numbers.length) {
        lines.push('');
        lines.push('| Number | In the answer | Source |');
        lines.push('| --- | --- | --- |');
        prov.numbers.forEach((n) => lines.push(`| ${n.literal}${n.unit ? ` ${n.unit}` : ''} | ${n.snippet.replace(/\|/g, '\\|')} | ${n.source ?? '**unsourced** (derived or invented)'} |`));
    }
    lines.push('');
    lines.push('## Safety');
    lines.push('');
    const safety = record.safety;
    lines.push(`- Unsafe actions: ${safety.unsafe_actions.length ? safety.unsafe_actions.map((u) => `step ${u.step}: ${u.what}`).join('; ') : '0'}`);
    lines.push(`- Unknown-tool calls: ${safety.unknown_tool_calls} · failed calls: ${safety.failed_calls} · questions to the user: ${safety.questions_asked} (all skipped by the simulated user)`);
    lines.push(`- Skills read: ${safety.skills_read.length ? safety.skills_read.join(', ') : 'none'}`);
    lines.push('');
    lines.push('## Planted truth (never shown to the model)');
    lines.push('');
    lines.push(record.truth.summary);
    lines.push('');
    lines.push(fence({ features: record.truth.features, tolerated: record.truth.tolerated, symmetry: record.truth.symmetry }));
    // The judgment subset (plan item 1.6 ii): which situation the case tests.
    if (record.truth.judgment) {
        lines.push('');
        lines.push(`A judgment case (plan item 1.6 ii, dev only): situation \`${record.truth.judgment.situation}\`.`);
    }
    if (record.truth.traps?.length) {
        lines.push('');
        lines.push('Traps in this case:');
        record.truth.traps.forEach((trap) => lines.push(`- ${trap}`));
    }
    lines.push('');
    return lines.join('\n');
};
