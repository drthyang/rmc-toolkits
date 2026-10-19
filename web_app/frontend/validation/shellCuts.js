// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// The two readings of whether a bond window cuts a coordination shell, side by
// side (architecture candidate #6). Both live in the run's one shell reading
// (src/llm/context/shells.js) and both are today's rules, unchanged:
//
//   asks     the Workbench's CN margin (window().cnMargin.cut): the window
//            overlaps the shell and holds fewer partner atoms than the whole
//            shell by at least max(0.5, 15 % of its CN) — clarify.js asks;
//   refuses  compute_bond_angles' counts (countsCut on the bond-angle
//            engine's coordinationSummary): more than half the central atoms
//            unbonded, or fewer than 60 % with the most common count.
//
// Which rule should decide a cut is the user's decision (the CN-margin rule
// for the tool is a recommendation: validation/findings/decisions/
// 2026-10-19-shell-catalogue.json). This module computes the evidence on a
// run's natural windows: for every closed shell, the triplet around the
// polyhedron's centre (the less abundant element of a mixed pair, as arm A
// reads a centre; the element itself for a like pair) in its whole shell, its
// first peak alone (two or more resolved peaks), and its top moved 0.2 Å
// either way. `validation/run.js shell-cuts` prints it for the dev cases and
// the Demo; the held-out test cases are never read.

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
 * Both readings on every natural window of one run: [{ run, triplet, label,
 * window, cn, shellCn, holds, asks, refuses, by, modal, modalShare,
 * unbondedShare }]; an empty window (a shell narrower than the 0.2 Å cut)
 * has `empty: true` and nothing else. `run` = { name, props, engine } (the
 * run module's CopilotRun: props.plotFiles, props.structure,
 * props.runSettings; engine(kind, params) the bond-angle engine).
 * `extraWindows` adds windows per "apex|end" ([label, [lo, hi]]).
 */
export const shellCutRows = async ({ name, props, engine }, { extraWindows = {} } = {}) => {
    const catalogue = shellCatalogue(props);
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
            const reading = catalogue.window(apex, end, window);
            const payload = await engine('triplets', { end1: end, apex, end2: end, r12Min: window[0], r12Max: window[1], binWidth: 1 });
            const summary = payload.coordinationSummary;
            const counts = countsCut(summary);
            rows.push({
                run: name,
                triplet,
                label,
                window: [...window],
                cn: round(reading.cn),
                shellCn: round(shell.cn),
                holds: reading.holds,
                asks: reading.cnMargin.cut,
                refuses: counts.cut,
                by: counts.by,
                modal: summary.modal,
                modalShare: round(summary.modalShare),
                unbondedShare: round(summary.unbondedShare)
            });
        }
    }
    return rows;
};

/** Where the readings agree and where they do not: { windows, both, neither, countsOnly, marginOnly, empty }. */
export const shellCutSummary = (rows) => {
    const read = rows.filter((row) => !row.empty);
    return {
        windows: read.length,
        both: read.filter((row) => row.asks && row.refuses),
        neither: read.filter((row) => !row.asks && !row.refuses),
        countsOnly: read.filter((row) => row.refuses && !row.asks),
        marginOnly: read.filter((row) => row.asks && !row.refuses),
        empty: rows.filter((row) => row.empty)
    };
};

const percent = (share) => `${(share * 100).toFixed(1)} %`;
const line = (row) => `| ${row.run} | ${row.triplet} | ${row.label} | ${row.window.join('–')} Å | ${row.cn} of ${row.shellCn} `
    + `| ${row.modal} at ${percent(row.modalShare)}${row.unbondedShare ? ` (${percent(row.unbondedShare)} unbonded)` : ''} |`;

/** The summary as Markdown: the totals, then every window the readings disagree on, and those both cut. */
export const renderShellCuts = (rows) => {
    const summary = shellCutSummary(rows);
    const table = (title, list) => (list.length
        ? [`**${title}** (${list.length})`, '', '| run | triplet | window | Å | CN of the shell\'s | modal count |', '|---|---|---|---|---|---|', ...list.map(line), '']
        : [`**${title}**: none`, '']);
    return [
        `${summary.windows} windows (${summary.empty.length} empty): both readings agree on ${summary.both.length + summary.neither.length} `
            + `(${summary.both.length} cut by both, ${summary.neither.length} by neither).`,
        '',
        ...table('compute_bond_angles refuses, the Workbench would not ask (counts only)', summary.countsOnly),
        ...table('The Workbench asks, compute_bond_angles would not refuse (CN margin only)', summary.marginOnly),
        ...table('Both', summary.both)
    ].join('\n');
};
