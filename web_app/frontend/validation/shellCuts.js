// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// The cut rule on a run's natural windows (architecture candidate #6). Since
// the user's approval (2026-10-19: validation/findings/decisions/
// 2026-10-19-shell-cut-rule.json) one rule decides whether a bond window cuts
// a coordination shell, the shell module's windowCut (src/llm/context/
// shells.js): where the window meets a shell, the CN margin — it overlaps the
// shell and holds fewer partner atoms than the whole shell by at least
// max(0.5, 15 % of its CN). For every window this module records
//
//   asks     the Workbench's question before the call (clarify.js beforeCall
//            on the tool's own windows);
//   refuses  compute_bond_angles' refusal (the tool as the run module builds
//            it, run on the window) — one rule, so the two agree;
//   before   the coordination counts' reading (countsCut on the bond-angle
//            engine's coordinationSummary: more than half the central atoms
//            unbonded, or fewer than 60 % with the most common count), the
//            tool's rule before the switch, kept as information: where it
//            differs from `refuses`, the switch changed the verdict.
//
// The natural windows: for every closed shell, the triplet around the
// polyhedron's centre (the less abundant element of a mixed pair, as arm A
// reads a centre; the element itself for a like pair) in its whole shell, its
// first peak alone (two or more resolved peaks), and its top moved 0.2 Å
// either way. `validation/run.js shell-cuts` prints it for the dev cases and
// the Demo; the held-out test cases are never read.

import { buildClarifier } from '../src/llm/copilot/clarify';
import { countsCut, shellCatalogue } from '../src/llm/context/shells';

const round = (value, digits = 3) => (Number.isFinite(value) ? Number(value.toPrecision(digits)) : null);

/** The natural triplet of each closed shell: [{ apex, end, shell }] (shell = the centre's own reading). */
export const naturalTriplets = (catalogue, elementCounts = {}) => (catalogue?.pairs || [])
    .filter((shell) => shell.window && Number.isFinite(shell.cn))
    .map((shell) => {
        const [apex, end] = shell.centre === shell.partner || (elementCounts[shell.centre] || 0) <= (elementCounts[shell.partner] || 0)
            ? [shell.centre, shell.partner]
            : [shell.partner, shell.centre];
        return { apex, end, shell: catalogue.shell(apex, end) };
    })
    .filter(({ shell }) => shell?.window && Number.isFinite(shell.cn));

/**
 * The cut rule on every natural window of one run: [{ run, triplet, label,
 * window, cn, shellCn, holds, asks, refuses, before, by, modal, modalShare,
 * unbondedShare }]; an empty window (a shell narrower than the 0.2 Å cut) has
 * `empty: true` and nothing else. `run` = { name, props, engine, tools } (the
 * run module's CopilotRun: props.plotFiles, props.structure,
 * props.runSettings; engine(kind, params) the bond-angle engine; tools its
 * engine tools, compute_bond_angles among them). `extraWindows` adds windows
 * per "apex|end" ([label, [lo, hi]]).
 */
export const shellCutRows = async ({ name, props, engine, tools }, { extraWindows = {} } = {}) => {
    const catalogue = shellCatalogue(props);
    const clarifier = buildClarifier(props);
    const tool = (tools || []).find((candidate) => candidate.name === 'compute_bond_angles');
    if (!tool) throw new Error(`${name}: no compute_bond_angles among the run's tools`);
    const rows = [];
    for (const { apex, end, shell } of naturalTriplets(catalogue, props.structure?.elementCounts)) {
        const windows = [
            ['whole shell', shell.window],
            ...(shell.firstPeak ? [['first peak', shell.firstPeak.window]] : []),
            ['cut 0.2 Å earlier', shell.cut.earlier.window],
            ['cut 0.2 Å later', shell.cut.later.window],
            ...(extraWindows[`${apex}|${end}`] || [])
        ];
        const triplet = `${end}–${apex}–${end}`;
        for (const [label, window] of windows) {
            if (!(window[1] > window[0])) {
                rows.push({ run: name, triplet, label, window: [...window], empty: true });
                continue;
            }
            const args = { end1: end, apex, end2: end, r12_min: window[0], r12_max: window[1] };
            const reading = catalogue.window(apex, end, window);
            const asks = clarifier.beforeCall(tool, args, { questionText: '' })?.kind === 'cut';
            const result = await tool.run(args);
            const summary = (await engine('triplets', { end1: end, apex, end2: end, r12Min: window[0], r12Max: window[1], binWidth: 1 })).coordinationSummary;
            const counts = countsCut(summary);
            rows.push({
                run: name,
                triplet,
                label,
                window: [...window],
                cn: round(reading.cn),
                shellCn: round(shell.cn),
                holds: reading.holds,
                asks,
                refuses: typeof result.data.refused === 'string',
                before: counts.cut,
                by: counts.by,
                modal: summary.modal,
                modalShare: round(summary.modalShare),
                unbondedShare: round(summary.unbondedShare)
            });
        }
    }
    return rows;
};

/**
 * One rule, and what it changed: { windows, refused, disagree (the Workbench
 * and the tool), changed, nowRefused (the CN margin cuts; the counts kept
 * them), noLongerRefused (the counts refused; the window holds the shell
 * within the margin), both (refused before and now), empty }.
 */
export const shellCutSummary = (rows) => {
    const read = rows.filter((row) => !row.empty);
    return {
        windows: read.length,
        refused: read.filter((row) => row.refuses),
        disagree: read.filter((row) => row.asks !== row.refuses),
        changed: read.filter((row) => row.refuses !== row.before),
        nowRefused: read.filter((row) => row.refuses && !row.before),
        noLongerRefused: read.filter((row) => !row.refuses && row.before),
        both: read.filter((row) => row.refuses && row.before),
        empty: rows.filter((row) => row.empty)
    };
};

const percent = (share) => `${(share * 100).toFixed(1)} %`;
const line = (row) => `| ${row.run} | ${row.triplet} | ${row.label} | ${row.window.join('–')} Å | ${row.cn} of ${row.shellCn} `
    + `| ${row.modal} at ${percent(row.modalShare)}${row.unbondedShare ? ` (${percent(row.unbondedShare)} unbonded)` : ''} |`;

/** The summary as Markdown: one rule's totals, then every window whose verdict the switch changed, and those refused by both. */
export const renderShellCuts = (rows) => {
    const summary = shellCutSummary(rows);
    const table = (title, list) => (list.length
        ? [`**${title}** (${list.length})`, '', '| run | triplet | window | Å | CN of the shell\'s | modal count |', '|---|---|---|---|---|---|', ...list.map(line), '']
        : [`**${title}**: none`, '']);
    return [
        `${summary.windows} windows (${summary.empty.length} empty): compute_bond_angles refuses ${summary.refused.length} and the Workbench `
            + `asks on the same windows (one rule, windowCut): they disagree on ${summary.disagree.length}. Against the counts rule the tool `
            + `refused on before 2026-10-19, ${summary.changed.length} change verdict: ${summary.noLongerRefused.length} no longer refused, `
            + `${summary.nowRefused.length} refused now (${summary.both.length} refused by both).`,
        '',
        ...table('The Workbench and compute_bond_angles disagree', summary.disagree),
        ...table('No longer refused (the counts refused; the window holds the shell within the CN margin)', summary.noLongerRefused),
        ...table('Refused now (the CN margin cuts; the counts kept them)', summary.nowRefused),
        ...table('Refused by both rules', summary.both)
    ].join('\n');
};
