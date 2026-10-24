// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

import { buildRunContext, seriesStats } from '../context/runContext';
import {
    BOUNDARY_REACH_A, CONSTRAINT_LIMIT_FRACTION, PEAK_SEARCH_MAX_A, averageShellCheck, averageStructureBonds, extractPeaks, shellCatalogue
} from '../context/shells';
import { MIN_TOTAL_DROP, TREND_MIN_ROWS, TREND_NOISE_FACTOR, chiSeries, judgeConvergence } from '../watchdog/heuristics';
import { oneOf } from './dataShape';
import { methodCheck } from './method';

// Copilot tools over the run data this module already receives as props
// (plotFiles, rValueFile, structure, symmetry, runSettings, pcaSites) — so
// they need nothing from the host. The engine tools (sites, displacement
// directions, bond angles, symmetry at a tolerance) come from the host as the
// `tools` prop. Each tool is offered only when its data exists for this run.
// Given the host's card catalog, each tool also names the cards it reads
// (`cards(args)`), which breathe while it runs: a card's `kind` (a plot kind,
// 'r_value', 'symmetry', 'model') and `file` tie it to the data. Each declares
// the shape of its `data` (`dataShape`, dataShape.js): the paths the model
// quotes by reference. Method advice in a description or a check is marked
// (`method`, method.js), so the validation's arm C reads none.

const roundSig = (value, digits = 3) => (
    Number.isFinite(value) && value !== 0 ? Number(value.toPrecision(digits)) : (Number.isFinite(value) ? value : null)
);

const SECTION_NOTES = {
    structure: 'cell lengths/angles of one unit cell, supercell, atom counts',
    configuration_optimization: 'move counters from the .rmc6f header',
    run_settings: 'the RMCProfile .dat control file: constraints, move sizes, fitted data',
    symmetry: 'detected space group, tolerance ladder, Wyckoff orbits with displacement',
    pca_displacements: 'per-site PCA of the displacement clouds (only after the PCA Ellipsoid page loaded them)',
    pair_correlations: 'partial g(r) first peaks vs the average-structure neighbour distance',
    datasets: 'per-dataset Rwp',
    convergence: 'ln(chi^2) history of the last .log column; the last-tenth trend of every chi^2 term and their total'
};

// Each section's shape as buildRunContext (context/runContext.js) builds it
// and read_run_data returns it, untrimmed.
const XYZ = ['number|null', 'number|null', 'number|null'];
const SECTION_SHAPES = {
    structure: {
        cell_A: XYZ,
        angles_deg: XYZ,
        supercell: ['number', 'number', 'number'],
        total_atoms: 'number',
        'elements?': { '*': 'number' }
    },
    configuration_optimization: {
        'moves_generated?': 'number',
        'moves_tried?': 'number',
        'moves_accepted?': 'number',
        'acceptance_ratio?': 'number|null',
        'accepted_moves_per_atom?': 'number|null',
        'accumulated_time_h?': 'number|null'
    },
    run_settings: {
        'title?': 'string',
        'material?': 'string',
        'phase?': 'string',
        'temperature?': 'string',
        // Labelled by element pair (atom) when the counts line up, else as listed.
        'min_distances_A?': oneOf({ '*-*': 'number' }, ['number']),
        'max_move_A?': oneOf({ '*': 'number' }, ['number']),
        'time_limit?': 'string',
        'save_period?': 'string',
        'weight_optimization?': 'boolean',
        'flags?': ['string'],
        'fitted_data?': [{ type: 'string', 'file?': 'string', 'fit?': 'string', 'fitted_in?': 'string' }],
        'independent_datasets?': 'number'
    },
    symmetry: {
        space_group: 'string',
        'point_group?': 'string',
        'n_ops?': 'number',
        'tolerance_A?': 'number',
        'max_residual_A?': 'number',
        'note?': 'string',
        'ladder?': [{ sg: 'string', 'pg?': 'string', holds_A: ['number|null', 'number|null'], n_ops: 'number' }],
        'point_group_facts?': ['string'],
        'noise_floor_A?': 'number',
        'noise_floor_note?': 'string',
        'sites?': [{
            element: 'string',
            multiplicity: 'number',
            'wyckoff?': 'string',
            'site_sym?': 'string',
            'frac?': XYZ,
            'mean_disp_A?': 'number',
            'max_disp_A?': 'number'
        }],
        'sites_omitted?': 'number'
    },
    pca_displacements: {
        note: 'string',
        sites: [{
            ref: 'number',
            element: 'string',
            U_iso_A2: 'number|null',
            'rms_axes_A?': XYZ,
            anisotropy: 'number|null',
            'non_gaussianity?': 'number',
            'degenerate?': 'boolean',
            'zero_spread?': 'boolean',
            'mixed?': 'boolean',
            'element_counts?': { '*': 'number' }
        }],
        'sites_omitted?': 'number'
    },
    pair_correlations: [{
        pair: 'string',
        'avg_structure_d_A?': 'number',
        'gr_peaks_A?': [{ r: 'number|null', fwhm: 'number|null' }],
        'first_shell?': {
            mean_d_A: 'number',
            to_A: 'number',
            'cn_*_per_*': 'number',
            'first_peak_to_A?': 'number|null',
            'first_peak_cn?': 'number|null',
            'gap_g_fraction?': 'number|null'
        },
        'cut_at_min_distance?': { at_A: 'number', g_fraction: 'number|null' }
    }],
    datasets: [{ kind: 'string', 'title?': 'string', 'rwp?': 'number', 'same_measurement_as?': 'string' }],
    convergence: {
        quantity: 'string',
        n_steps: 'number',
        first: 'number|null',
        last: 'number|null',
        min: 'number|null',
        max: 'number|null',
        watchdog_status: 'string',
        history: ['number|null'],
        // Each series' fitted ln(chi^2) change over the last tenth, the
        // scatter about the fit and its trend (watchdog/heuristics.js); the
        // total only with the log's terms.
        'last_10pct?': {
            rows: 'number',
            'total?': { change: 'number|null', noise_sd: 'number|null', trend: 'string' },
            terms: { '*': { change: 'number|null', noise_sd: 'number|null', trend: 'string' } },
            'zero_terms?': 'number'
        },
        'non_finite_steps?': 'number',
        'non_finite_note?': 'string',
        'column?': 'string',
        'final_chi_squared?': 'number'
    }
};

const pairKey = (label) => {
    const match = /^([A-Z][a-z]?)\s*-\s*([A-Z][a-z]?)$/.exec((label || '').trim());
    return match ? [match[1], match[2]].sort().join('-') : null;
};

// The card kinds a read_run_data section is about.
const SECTION_CARD_KINDS = {
    structure: ['model'],
    configuration_optimization: ['model'],
    symmetry: ['symmetry'],
    pair_correlations: ['pdf_partials'],
    convergence: ['r_value']
};

const readRunDataTool = (context, cardIds) => {
    const sections = Object.keys(SECTION_NOTES).filter((key) => context[key] != null);
    if (!sections.length) return null;
    return {
        name: 'read_run_data',
        description: 'Read one section of the run summary in full (the context you were given is trimmed '
            + 'to fit). Sections: '
            + sections.map((key) => `${key} (${SECTION_NOTES[key]})`).join('; ') + '.',
        parameters: {
            type: 'object',
            properties: { section: { type: 'string', enum: sections, description: 'Which section to read.' } },
            required: ['section']
        },
        // The section's name, and the section under that name (one of these).
        dataShape: {
            section: 'string',
            ...Object.fromEntries(sections.map((key) => [`${key}?`, SECTION_SHAPES[key]]))
        },
        label: ({ section }) => `Read ${String(section || '').replace(/_/g, ' ')}`,
        cards: ({ section }) => cardIds((card) => (SECTION_CARD_KINDS[section] || []).includes(card.kind)),
        run: async ({ section }) => ({ data: { section, [section]: context[section] }, checks: [] })
    };
};

// "a", "a and b", "a, b and c".
const listOf = (names) => (names.length > 1 ? `${names.slice(0, -1).join(', ')} and ${names.at(-1)}` : names[0]);

// One series' last-tenth verdict as the tool gives it (heuristics.js seriesTrend).
const trendRow = (entry) => ({
    change_last_10pct: roundSig(entry.change, 2),
    noise_sd: roundSig(entry.noiseSd, 2),
    trend: entry.trend
});
const TREND_ROW = { change_last_10pct: 'number|null', noise_sd: 'number|null', trend: 'string' };

// Why a run is diverging, from the verdict (heuristics.js judgeConvergence).
const divergingBecause = (judgement, nonFiniteLatest) => {
    if (nonFiniteLatest) return 'the latest log row is non-finite (NaN/Inf): the run blew up';
    if (judgement.total?.trend === 'rising') return `the total rises by more than ${TREND_NOISE_FACTOR} × its noise over the last tenth`;
    const rising = judgement.rising.filter((name) => name !== 'total');
    return `${listOf(rising)} ${rising.length > 1 ? 'rise' : 'rises'} by more than ${TREND_NOISE_FACTOR} × ${rising.length > 1 ? 'their' : 'its'} `
        + 'noise over the last tenth while nothing falls';
};

const convergenceTool = (rValueFile, cardIds) => {
    const history = rValueFile?.plotData?.series?.[0]?.y;
    if (!Array.isArray(history) || !history.length) return null;
    const column = rValueFile?.plotData?.chiColumn || null;
    return {
        name: 'assess_convergence',
        description: 'Classify the convergence of every chi^2 term of the .log and of their total (ln chi^2 per log row; '
            + `the dashboard charts only the last column), with the dashboard watchdog's rule: a series is still falling when `
            + `the least-squares line through its last tenth (at least ${TREND_MIN_ROWS} rows, never back past the last restart) `
            + `drops by more than ${TREND_NOISE_FACTOR} × the scatter of its points about that line; the run is converged only when `
            + 'no term and not the total still falls (or rises), stalled when it is flat without any of them ever dropping by 10 %. '
            + 'Status: improving | converged | stalled | diverging.',
        parameters: { type: 'object', properties: {} },
        // total only with the log's chi^2 terms (a log whose line 1 names its
        // columns); without them the charted column is the one term.
        dataShape: {
            quantity: 'string',
            status: 'string',
            n_steps: 'number',
            window_rows: 'number',
            'total?': { first: 'number|null', last: 'number|null', min: 'number|null', drop: 'number|null', ...TREND_ROW },
            terms: [{ term: 'string', last: 'number|null', ...TREND_ROW }],
            'zero_terms?': ['string'],
            chart_column: 'string|null',
            'non_finite_steps?': 'number'
        },
        label: () => 'Convergence check',
        page: 'dashboard',
        cards: () => cardIds((card) => card.kind === 'r_value'),
        run: async () => {
            const stats = seriesStats(history);
            const series = chiSeries(rValueFile.plotData);
            const judgement = judgeConvergence(series);
            const { status } = judgement;
            const data = {
                quantity: 'ln(chi^2) per .log row, natural log; lower is better',
                status,
                n_steps: stats.nSteps,
                window_rows: judgement.windowRows
            };
            if (judgement.total) {
                const { first, last, min, drop } = judgement.total;
                data.total = { first: roundSig(first), last: roundSig(last), min: roundSig(min), drop: roundSig(drop), ...trendRow(judgement.total) };
            }
            data.terms = judgement.terms.map((term) => ({ term: term.name, last: roundSig(term.last), ...trendRow(term) }));
            if (judgement.zeroTerms.length) data.zero_terms = judgement.zeroTerms;
            data.chart_column = column;
            if (stats.nonFiniteSteps) data.non_finite_steps = stats.nonFiniteSteps;
            const noisy = `${TREND_NOISE_FACTOR} × the noise`;
            const checks = [{
                level: 'note',
                text: judgement.total
                    ? `Judged on every chi^2 term of the .log (${judgement.terms.length} non-zero`
                        + `${judgement.zeroTerms.length ? `, ${judgement.zeroTerms.length} zero throughout` : ''}) and their total, over the last `
                        + `tenth of the run (${judgement.windowRows} rows): converged only when none falls beyond ${noisy}.`
                    : `Judged on ${column ? `the ${column} term` : 'the last .log column'} only: the log names no chi^2 columns, so no other `
                        + 'term and no total is known.'
            }];
            const nonFiniteLatest = !Number.isFinite(history[history.length - 1])
                || (Array.isArray(series.total) && !Number.isFinite(series.total[series.total.length - 1]))
                || series.terms.some((term) => !Number.isFinite(term.values[term.values.length - 1]));
            if (status === 'unknown') {
                checks.push({ level: 'fail', text: 'Fewer than three finite log rows in the last tenth — no trend to judge.' });
            } else if (status === 'converged') {
                const best = [...judgement.terms, ...(judgement.total ? [{ ...judgement.total, name: 'the total' }] : [])]
                    .reduce((top, entry) => ((entry.drop ?? -Infinity) > (top?.drop ?? -Infinity) ? entry : top), null);
                checks.push({
                    level: 'ok',
                    text: `Flat: no term${judgement.total ? ' and not the total' : ''} falls or rises beyond ${noisy} over the last tenth, `
                        + `after ${best.name} dropped by ${roundSig(best.drop)} in ln(chi^2).`
                });
            } else if (status === 'improving') {
                const falling = judgement.falling.map((name) => (name === 'total' ? 'the total' : name));
                checks.push({
                    level: 'note',
                    text: `Still improving: ${listOf(falling)} ${falling.length > 1 ? 'fall' : 'falls'} by more than ${TREND_NOISE_FACTOR} × `
                        + `${falling.length > 1 ? 'their' : 'its'} noise over the last tenth.`
                });
            } else if (status === 'stalled') {
                checks.push({
                    level: 'warn',
                    text: `Stalled: flat over the last tenth, but no term${judgement.total ? ' and not the total' : ''} ever dropped by `
                        + `${MIN_TOTAL_DROP} in ln(chi^2) (10 %).`
                });
            } else {
                checks.push({ level: 'warn', text: `Diverging: ${divergingBecause(judgement, nonFiniteLatest)}.` });
            }
            if (stats.nonFiniteSteps) checks.push({ level: 'warn', text: `${stats.nonFiniteSteps} log rows are non-finite (NaN/Inf).` });
            if (stats.nSteps < 5) checks.push({ level: 'warn', text: `Only ${stats.nSteps} log rows.` });
            // For the answer check: never "converged" against this status.
            return { data, checks, facts: { convergence: status } };
        }
    };
};

// The shell block, the first peak, the checks, the facts and the constraint
// come from the run's one shell reading (context/shells.js); this tool only
// shows it. Its peaks_A list is the caller's display, below r_max.
const grPeaksTool = (run, cardIds) => {
    const partials = (run?.plotFiles || []).find(
        (file) => file.plotKind === 'pdf_partials' && file.plotData?.series?.length
    );
    if (!partials) return null;
    const series = partials.plotData.series.filter((entry) => pairKey(entry.label));
    if (!series.length) return null;
    const labels = series.map((entry) => entry.label.trim());
    const catalogue = shellCatalogue(run);
    return {
        name: 'find_gr_peaks',
        description: 'Analyse one partial g(r) (from the run\'s partial-PDF file): its peaks below r_max '
            + '(position, approximate FWHM) and its first coordination shell, which does not depend on r_max — the '
            + 'window from onset to the first deep minimum (it holds '
            + 'the whole shell, split shells included), mean distance and spread from r²g(r), coordination numbers in '
            + 'both directions with their sensitivity to the cut, whether the shell is separated from the next, '
            + 'for a split shell joined across a deep trough its parts (CN, share, mean distance) and the long − short '
            + 'bond difference, the average structure\'s bond lengths for the pair against the window, and whether the pair\'s '
            + 'hard minimum-distance constraint cuts the distribution.',
        // The method advice (method.js), inside the window's parenthesis: sent with method 'full'.
        method: [{ after: 'to the first deep minimum (', text: 'use first_shell.window_A for compute_bond_angles: ' }],
        parameters: {
            type: 'object',
            properties: {
                pair: { type: 'string', enum: labels, description: 'Element pair, as labelled in the partials file.' },
                r_max: { type: 'number', minimum: 1, maximum: 20, description: 'List the peaks up to this r (Å). Default 6.' },
                max_peaks: { type: 'integer', minimum: 1, maximum: 8, description: 'Default 3.' }
            },
            required: ['pair']
        },
        // No peak below r_max and no shell: pair, r_max_A and peaks_A [] only. The cn keys name the
        // count each way, <partner>_per_<centre> (Se_per_Ta). The first peak
        // alone (first_peak_only) only for a shell of two or more resolved
        // peaks; its cn only with the cell. The split (plan item 4.2) only for
        // a shell joined across a deep trough (shells.js): its parts, short
        // to long, each one's cn in the direction the join decided (one key),
        // and long − short with its sampling error.
        dataShape: {
            pair: 'string',
            r_max_A: 'number',
            peaks_A: [{ r: 'number|null', fwhm: 'number|null' }],
            'first_shell?': {
                window_A: ['number|null', 'number|null'],
                peaks_in_shell_A: ['number|null'],
                'mean_d_A?': 'number|null',
                'sd_d_A?': 'number|null',
                'cn?': { '*_per_*': 'number|null' },
                'cn_if_cut_0p2A_earlier_later?': ['number|null', 'number|null'],
                trough_g_fraction: 'number|null',
                'first_peak_only?': {
                    window_A: ['number|null', 'number|null'],
                    'cn?': { '*_per_*': 'number|null' },
                    trough_g_fraction: 'number|null'
                },
                'split?': {
                    parts: [{
                        window_A: ['number|null', 'number|null'],
                        mean_d_A: 'number|null',
                        sd_d_A: 'number|null',
                        cn: { '*_per_*': 'number|null' },
                        share: 'number|null'
                    }],
                    long_minus_short_A: 'number|null',
                    long_minus_short_se_A: 'number|null'
                },
                // The average structure's first shell for the pair (shortest,
                // longest bond between site means); only with its basis.
                'average_structure_d_A?': ['number|null', 'number|null']
            },
            'min_distance_constraint_A?': 'number',
            'g_at_constraint_fraction?': 'number|null'
        },
        label: ({ pair }) => `g(r) peaks ${pair || ''}`.trim(),
        page: 'dashboard',
        cards: () => cardIds((card) => Boolean(partials.name) && card.file === partials.name),
        run: async ({ pair, r_max: rMax = 6, max_peaks: maxPeaks = 3 }) => {
            const entry = series.find((item) => item.label.trim() === pair);
            const peaks = extractPeaks(entry.x, entry.y, { rMax, maxPeaks });
            const data = { pair, r_max_A: rMax, peaks_A: peaks.map((peak) => ({ r: roundSig(peak.r, 4), fwhm: roundSig(peak.fwhm, 2) })) };
            const checks = [];
            const [a, b] = pair.split('-').map((element) => element.trim());
            // The first coordination shell, as the run reads it: peaks below
            // 8 Å, the closing minimum to max(8, r₁ + 3) Å, whatever r_max is.
            const shell = catalogue?.shell(a, b) ?? null;
            if (!peaks.length && !shell) {
                checks.push({ level: 'fail', text: `No peak above 15 % of the maximum below ${rMax} Å.` });
                return { data, checks };
            }
            if (!peaks.length) {
                checks.push({
                    level: 'note',
                    text: `No peak above 15 % of the maximum below ${rMax} Å (r_max); the first shell, read from the peaks below `
                        + `${PEAK_SEARCH_MAX_A} Å, lies beyond it.`
                });
            }
            const firstR = peaks.length ? peaks[0].r : shell.peaks[0].r;
            const perA = `${b}_per_${a}`;
            if (shell?.window) {
                const block = {
                    window_A: [...shell.window],
                    peaks_in_shell_A: shell.peaksInShell.map((r) => roundSig(r, 4))
                };
                if (Number.isFinite(shell.meanR)) {
                    block.mean_d_A = roundSig(shell.meanR, 4);
                    block.sd_d_A = roundSig(shell.sdR, 2);
                }
                if (Number.isFinite(shell.cn)) {
                    block.cn = { [perA]: roundSig(shell.cn) };
                    if (a !== b && Number.isFinite(shell.reverseCn)) block.cn[`${a}_per_${b}`] = roundSig(shell.reverseCn);
                    block.cn_if_cut_0p2A_earlier_later = [shell.cut.earlier.cn, shell.cut.later.cn].map((value) => roundSig(value));
                }
                block.trough_g_fraction = roundSig(shell.troughRatio, 2);
                // Two or more resolved peaks: the first alone, to the minimum
                // between the first two.
                if (shell.firstPeak) {
                    block.first_peak_only = {
                        window_A: [...shell.firstPeak.window],
                        ...(Number.isFinite(shell.firstPeak.cn) ? { cn: { [perA]: roundSig(shell.firstPeak.cn) } } : {}),
                        trough_g_fraction: roundSig(shell.firstPeak.troughRatio, 2)
                    };
                }
                // A shell joined across a deep trough: its split by number
                // (shells.js, plan item 4.2), each part's count in the
                // direction the join decided — the one whose sum is ideal.
                const { split } = shell;
                const forward = split?.direction !== 'reverse';
                const [counted, around] = forward ? [b, a] : [a, b];
                const splitKey = `${counted}_per_${around}`;
                if (split) {
                    block.split = {
                        parts: split.parts.map((part) => ({
                            window_A: [...part.window],
                            mean_d_A: roundSig(part.meanR, 4),
                            sd_d_A: roundSig(part.sdR, 2),
                            cn: { [splitKey]: roundSig(forward ? part.cn : part.reverseCn) },
                            share: roundSig(part.share, 3)
                        })),
                        long_minus_short_A: roundSig(split.longMinusShort, 3),
                        long_minus_short_se_A: roundSig(split.longMinusShortSe, 2)
                    };
                }
                // The average structure's first shell for the pair, between site
                // means (shells.js averageStructureBonds; the user's rule,
                // 2026-10-24): the window is cross-checked against it.
                const average = averageShellCheck(shell, averageStructureBonds(run?.structure, a, b));
                if (average) block.average_structure_d_A = average.bonds.map((r) => roundSig(r, 3));
                data.first_shell = block;
                if (average?.disagrees) {
                    // What to do about it is method (method.js): arm C's method 'off' reads the reading alone.
                    checks.push(methodCheck('warn',
                        `The average structure's ${a}–${b} bonds (${block.average_structure_d_A.join('–')} Å between the site means, `
                            + `first_shell.average_structure_d_A) lie up to ${roundSig(average.outsideA, 2)} Å outside first_shell.window_A `
                            + `(${shell.window[0]}–${shell.window[1]} Å), more than the shell's width (its first peak's FWHM, `
                            + `${roundSig(average.widthA, 2)} Å): the configuration's first shell and the average structure's disagree`,
                        { method: '; check the window before using it' },
                        '.'));
                }
                const inShell = shell.peaksInShell.map((r) => roundSig(r, 4));
                if (block.split && shell.firstPeak) {
                    const { parts, long_minus_short_A: longMinusShort, long_minus_short_se_A: se } = block.split;
                    const counts = parts.map((part, k) => `${part.cn[splitKey]} ${k ? '' : `${counted} per ${around} `}at ${part.mean_d_A} Å`);
                    const bonds = parts.length === 2 ? 'the short and long bonds' : `the ${parts.length} bond lengths, short to long,`;
                    // Which window to use is method (method.js): arm C's method 'off' sends the reading alone.
                    checks.push(methodCheck('note',
                        `One split shell across a deep trough: ${counts.join(' + ')} (first_shell.split.parts), `
                            + `${roundSig(forward ? shell.cn : shell.reverseCn)} together, an ideal polyhedral count: ${bonds} of one distorted `
                            + `polyhedron (e.g. Jahn–Teller). Long − short ${longMinusShort}${se != null ? ` ± ${se}` : ''} Å `
                            + `(first_shell.split: long_minus_short_A${se != null ? ' ± long_minus_short_se_A, the configuration\'s sampling error' : ''}); `
                            + `${parts.map((part) => `${Math.round(part.share * 100)} %`).join(' / ')} of the shell's count. `
                            + 'first_shell.window_A holds every part; the first peak alone (first_peak_only, to '
                            + `${roundSig(shell.firstPeak.window[1], 3)} Å) cuts the polyhedron`,
                        { method: '; use it only if the user chooses it' },
                        '.'));
                }
                if (!shell.separated) {
                    // What to ask about the gap is method (method.js): arm C's method 'off' sends the reading alone.
                    const head = `The first shell is not separated from the next: at its boundary (${roundSig(shell.window[1], 3)} Å) `
                        + `g(r) is still ${Math.round(shell.troughRatio * 100)} % of the peak, so ${a}–${b} distances fill the gap `
                        + 'between shells';
                    checks.push(Number.isFinite(shell.cn)
                        ? methodCheck('warn',
                            `${head} and the coordination number depends on the cut (${roundSig(shell.cut.earlier.cn)}–${roundSig(shell.cut.later.cn)} `
                                + 'for a cut 0.2 Å earlier or later).',
                            { method: ` Ask whether the crystal has ${a}–${b} distances there.` })
                        : { level: 'warn', text: `${head}.` });
                }
                if (shell.peaksInShell.length > 1 && !shell.joined) {
                    const alone = Number.isFinite(shell.firstPeak?.cn) && Number.isFinite(shell.cn)
                        ? ` The first peak alone (first_peak_only, to ${roundSig(shell.firstPeak.window[1], 3)} Å) holds `
                            + `${roundSig(shell.firstPeak.cn)} ${b} per ${a} of the shell's ${roundSig(shell.cn)}.`
                        : '';
                    checks.push(methodCheck('note',
                        `The first shell holds ${shell.peaksInShell.length} peaks (${inShell.join(', ')} Å): `
                            + 'a split shell, e.g. two bond lengths of one distorted polyhedron.',
                        { method: ' Use first_shell.window_A, not the first peak alone.' },
                        alone));
                }
            } else if (shell) {
                const reach = Math.max(PEAK_SEARCH_MAX_A, shell.peaks[0].r + BOUNDARY_REACH_A);
                checks.push({ level: 'note', text: `No minimum closes the first shell below ${roundSig(reach, 3)} Å.` });
            } else {
                checks.push({ level: 'note', text: `No first shell is read: no peak below ${PEAK_SEARCH_MAX_A} Å, where the shells are read.` });
            }

            if (shell?.constraint) {
                const { limitA: limit, gFraction, limited } = shell.constraint;
                data.min_distance_constraint_A = limit;
                // g just above the constraint (its first 2.5 grid steps),
                // relative to the first peak: a cut-off tail piles up there.
                // Above CONSTRAINT_LIMIT_FRACTION (5 %) the peak is
                // constraint-limited (the user's rule, 2026-10-24). What to do
                // about it is method (method.js): arm C's method 'off' reads
                // the reading alone.
                if (limited) {
                    data.g_at_constraint_fraction = roundSig(gFraction, 2);
                    checks.push(methodCheck('warn',
                        `The ${limit} Å closest-approach constraint cuts this distribution: g(r) just above it is `
                            + `${Math.round(gFraction * 100)} % of the first peak (above ${Math.round(CONSTRAINT_LIMIT_FRACTION * 100)} %), `
                            + `so ${a}–${b} pairs pile up against the constraint and the constraint, not the data, shapes the low-r side `
                            + 'of the peak — constraint-limited',
                        { method: ': interpret nothing on its low-r side; consider rerunning with the constraint loosened' },
                        '.'));
                }
            }
            if (!checks.some((check) => check.level === 'warn')) {
                checks.push({
                    level: 'ok',
                    text: `First peak at ${roundSig(firstR, 4)} Å`
                        + (data.first_shell?.mean_d_A ? `; first shell to ${data.first_shell.window_A[1]} Å, mean ${data.first_shell.mean_d_A} Å` : '')
                        + (data.first_shell?.cn ? `, CN ${Object.entries(data.first_shell.cn).map(([k, v]) => `${v} ${k.replace(/_/g, ' ')}`).join(', ')}` : '')
                        + '.'
                });
            }
            // For the coordination answer check: the whole first shell's count,
            // each way (B around A, and A around B for a mixed pair).
            if (shell?.window && Number.isFinite(shell.cn)) {
                const window = data.first_shell.window_A;
                const counts = [{ centre: a, partner: b, window, cn: roundSig(shell.cn), whole_shell: true, confirmed: false }];
                const reverse = data.first_shell.cn[`${a}_per_${b}`];
                if (a !== b && Number.isFinite(reverse)) counts.push({ centre: b, partner: a, window, cn: reverse, whole_shell: true, confirmed: false });
                return { data, checks, facts: { coordination: counts } };
            }
            return { data, checks };
        }
    };
};

// Fit files with a model and a measured column, by plot kind.
const FIT_KINDS = {
    xpdf: 'X-ray G(r)',
    npdf: 'neutron G(r)',
    xray_sq: 'X-ray F(Q)',
    neutron_sq: 'neutron S(Q)',
    bragg: 'Bragg profile'
};
// RMCProfile writes (x, calculated, experimental); a label naming a role
// overrides the order. Duplicated from browserData.rwpColumns on purpose:
// this module may not import host code.
const EXPERIMENTAL_LABEL = /exp|obs/i;
const CALCULATED_LABEL = /calc|rmc|fit/i;
const fitSeriesRoles = (series) => {
    const experimental = series.findIndex((entry) => EXPERIMENTAL_LABEL.test(entry.label || ''));
    const calculated = series.findIndex((entry) => !EXPERIMENTAL_LABEL.test(entry.label || '') && CALCULATED_LABEL.test(entry.label || ''));
    return experimental >= 0 && calculated >= 0 ? [series[calculated], series[experimental]] : [series[0], series[1]];
};
// One range holding at least this share of the squared residual is where
// the misfit lives.
export const MISFIT_CONCENTRATION = 0.4;
// A range whose own Rw is this many times the dataset's is poorly fitted for
// its signal, whatever its share of the squared residual.
export const RELATIVE_MISFIT_FACTOR = 3;
// A residual feature is above noise when it is at least this many σ(r).
export const ABOVE_NOISE_SIGMA = 3;

// σ(r) of an X-ray G(r) (RMCProfile's transform of F(Q)): the file's own σ
// column, else the one propagated from F(Q)'s σ (plotData.uncertainty, the
// host's ft_uncertainty payload; docs/algorithms/run-dashboard.md, "The G(r)
// error ribbon"). Null for any other dataset, or a G(r) with no word on σ;
// { sigma: null, why } when its σ is known to be unavailable.
const fitSigma = (file, exp) => {
    if (file.plotKind !== 'xpdf') return null;
    if (Array.isArray(exp.sigma) && exp.sigma.length === exp.x.length) {
        return { sigma: 'file', from: file.name || null, values: exp.sigma, uncertainty: null };
    }
    const uncertainty = file.plotData?.uncertainty;
    if (uncertainty?.status === 'propagated' && Array.isArray(uncertainty.sigma) && uncertainty.sigma.length === exp.x.length) {
        return { sigma: 'propagated', from: uncertainty.sources?.input || null, values: uncertainty.sigma, uncertainty };
    }
    return uncertainty?.reason ? { sigma: null, why: uncertainty.reason } : null;
};

const fitTool = (plotFiles, cardIds) => {
    const files = (plotFiles || []).filter((file) => FIT_KINDS[file.plotKind] && file.plotData?.series?.length >= 2);
    if (!files.length) return null;
    const keyed = new Map();
    files.forEach((file) => {
        let key = file.plotKind;
        for (let n = 2; keyed.has(key); n += 1) key = `${file.plotKind}_${n}`;
        keyed.set(key, file);
    });
    return {
        name: 'compare_fit',
        description: 'Where the model misses the measured data, for one fitted dataset: the R-factor (as the '
            + 'dashboard\'s Rwp) over consecutive ranges of r or Q, the share of the squared residual in each range, '
            + 'and the largest residual features (position, model above or below the data, size relative to the '
            + 'data\'s rms). For an X-ray G(r) with σ(r) (propagated from F(Q)\'s σ), each feature\'s size in σ and '
            + 'whether it is above noise; its r points are correlated.',
        // The method advice (method.js): sent with method 'full', after the description.
        method: 'Use it for "where does the fit fail" and "is the misfit at short or long range"; run it once '
            + 'per dataset (a G(r) and its F(Q) are separate entries).',
        parameters: {
            type: 'object',
            properties: {
                dataset: {
                    type: 'string',
                    enum: [...keyed.keys()],
                    description: [...keyed.entries()].map(([key, file]) => `${key}: ${FIT_KINDS[file.plotKind]}`).join('; ')
                },
                n_ranges: { type: 'integer', minimum: 2, maximum: 12, description: 'Equal-width ranges over the data. Default 6.' }
            },
            required: ['dataset']
        },
        // Fewer than 10 points: the first five keys only. A range and a
        // residual are named by the axis: r_A / r in real space, Q_A / Q in
        // reciprocal space.
        dataShape: {
            dataset: 'string',
            kind: 'string',
            axis: 'string',
            model_series: 'string',
            data_series: 'string',
            'rw?': 'number|null',
            'range?': ['number|null', 'number|null'],
            'by_range?': [oneOf(
                { r_A: ['number|null', 'number|null'], rw: 'number|null', share_of_misfit: 'number|null' },
                { Q_A: ['number|null', 'number|null'], rw: 'number|null', share_of_misfit: 'number|null' }
            )],
            'largest_residuals?': [oneOf(
                {
                    r: 'number|null', model_minus_data: 'number|null', relative_to_data_rms: 'number|null',
                    'relative_to_sigma?': 'number|null', 'above_noise?': 'boolean'
                },
                { Q: 'number|null', model_minus_data: 'number|null', relative_to_data_rms: 'number|null' }
            )],
            // An X-ray G(r) only: its σ(r) and what the residual is against it,
            // or why it has none.
            'noise?': oneOf(
                {
                    sigma: 'string',
                    from: 'string|null',
                    rms_residual_over_sigma: 'number|null',
                    above_noise_at: 'number',
                    'reproduction?': { residual: 'number|null', tolerance: 'number|null' },
                    'correlation?': { neighbour: 'number|null', effective_points: 'number|null', points: 'number' }
                },
                { sigma: 'null', why: 'string' }
            )
        },
        label: ({ dataset }) => `Fit residual ${dataset || ''}`.trim(),
        page: 'dashboard',
        cards: ({ dataset }) => cardIds((card) => card.file && card.file === keyed.get(dataset)?.name),
        run: async ({ dataset, n_ranges: nRanges = 6 }) => {
            const file = keyed.get(dataset);
            const [calc, exp] = fitSeriesRoles(file.plotData.series);
            const points = [];
            exp.x.forEach((x, i) => {
                const observed = exp.y[i];
                const model = calc.y?.[i];
                if (Number.isFinite(x) && Number.isFinite(observed) && Number.isFinite(model)) points.push([x, observed, model, i]);
            });
            const axis = /sq|bragg/.test(file.plotKind) ? 'Q' : 'r';
            const data = { dataset, kind: FIT_KINDS[file.plotKind], axis, model_series: calc.label, data_series: exp.label };
            const checks = [];
            if (points.length < 10) {
                checks.push({ level: 'fail', text: `Only ${points.length} points with both a model and a measured value.` });
                return { data, checks };
            }
            const sumSq = (rows, f) => rows.reduce((sum, row) => sum + f(row) ** 2, 0);
            const totalResidual = sumSq(points, ([, observed, model]) => model - observed);
            const totalData = sumSq(points, ([, observed]) => observed);
            data.rw = roundSig(Math.sqrt(totalResidual / totalData), 3);
            data.range = [roundSig(points[0][0], 4), roundSig(points.at(-1)[0], 4)];
            const lo = points[0][0];
            const width = (points.at(-1)[0] - lo) / nRanges;
            data.by_range = Array.from({ length: nRanges }, (_, k) => {
                const from = lo + k * width;
                const to = k === nRanges - 1 ? Infinity : from + width;
                const rows = points.filter(([x]) => x >= from && x < to);
                const residual = sumSq(rows, ([, observed, model]) => model - observed);
                const signal = sumSq(rows, ([, observed]) => observed);
                return {
                    [`${axis}_A`]: [roundSig(from, 3), roundSig(Math.min(lo + (k + 1) * width, points.at(-1)[0]), 3)],
                    rw: signal > 0 ? roundSig(Math.sqrt(residual / signal), 2) : null,
                    share_of_misfit: roundSig(totalResidual > 0 ? residual / totalResidual : 0, 2)
                };
            });
            // The largest residual features: local extrema of a 5-point mean of
            // (model − data), relative to the data's rms over the whole range.
            const residuals = points.map(([, observed, model]) => model - observed);
            const smooth = residuals.map((_, i) => {
                let sum = 0;
                let n = 0;
                for (let k = -2; k <= 2; k++) {
                    if (i + k >= 0 && i + k < residuals.length) { sum += residuals[i + k]; n += 1; }
                }
                return sum / n;
            });
            const dataRms = Math.sqrt(totalData / points.length);
            const extrema = [];
            for (let i = 1; i < smooth.length - 1; i++) {
                const a = Math.abs(smooth[i]);
                if (a >= Math.abs(smooth[i - 1]) && a > Math.abs(smooth[i + 1])) extrema.push(i);
            }
            // An X-ray G(r)'s σ(r): a residual feature is above noise only
            // against it (≥ ABOVE_NOISE_SIGMA σ at its position; the 5-point
            // mean of correlated points has no smaller σ than one point).
            const sigma = fitSigma(file, exp);
            const sigmaAt = (i) => {
                const value = sigma?.values?.[points[i][3]];
                return Number.isFinite(value) && value > 0 ? value : null;
            };
            data.largest_residuals = extrema
                .sort((a, b) => Math.abs(smooth[b]) - Math.abs(smooth[a]))
                .slice(0, 5)
                .map((i) => {
                    const entry = {
                        [axis]: roundSig(points[i][0], 4),
                        model_minus_data: roundSig(smooth[i], 2),
                        relative_to_data_rms: roundSig(smooth[i] / dataRms, 2)
                    };
                    if (sigma?.values) {
                        const s = sigmaAt(i);
                        entry.relative_to_sigma = s === null ? null : roundSig(smooth[i] / s, 2);
                        entry.above_noise = s !== null && Math.abs(smooth[i]) >= ABOVE_NOISE_SIGMA * s;
                    }
                    return entry;
                })
                .sort((a, b) => a[axis] - b[axis]);
            if (sigma?.values) {
                const ratios = points.map((_, i) => {
                    const s = sigmaAt(i);
                    return s === null ? null : residuals[i] / s;
                }).filter((value) => value !== null);
                data.noise = {
                    sigma: sigma.sigma,
                    from: sigma.from,
                    rms_residual_over_sigma: ratios.length
                        ? roundSig(Math.sqrt(ratios.reduce((sum, value) => sum + value * value, 0) / ratios.length), 3)
                        : null,
                    above_noise_at: ABOVE_NOISE_SIGMA
                };
                const { uncertainty } = sigma;
                if (uncertainty) {
                    data.noise.reproduction = {
                        residual: roundSig(uncertainty.reproduction?.gr?.residual, 2),
                        tolerance: roundSig(uncertainty.reproduction?.tolerance, 2)
                    };
                    const { correlation } = uncertainty;
                    if (correlation) {
                        data.noise.correlation = {
                            neighbour: roundSig(correlation.neighbour, 2),
                            effective_points: Number.isFinite(correlation.effectivePoints) ? Math.round(correlation.effectivePoints) : null,
                            points: correlation.points
                        };
                    }
                }
            } else if (sigma?.why) {
                data.noise = { sigma: null, why: sigma.why };
            }
            const worst = data.by_range.reduce((best, entry, k) => (entry.share_of_misfit > data.by_range[best].share_of_misfit ? k : best), 0);
            const span = data.by_range[worst][`${axis}_A`];
            const unit = axis === 'r' ? ' Å' : ' Å⁻¹';
            // The finding leads the step card (its headline is the ok check).
            const where = data.by_range[worst].share_of_misfit >= MISFIT_CONCENTRATION
                ? `the misfit is concentrated in ${axis} = ${span[0]}–${span[1]}${unit}: `
                    + `${Math.round(data.by_range[worst].share_of_misfit * 100)} % of the squared residual (Rw ${data.by_range[worst].rw} there).`
                : `the misfit is spread over the range: no ${axis} range holds ${Math.round(MISFIT_CONCENTRATION * 100)} % of it.`;
            // A range where the misfit is large relative to its own (small)
            // signal: the share hides it — high-Q F(Q) carries the sharpest
            // local detail, and its oscillations are weak.
            const poorest = data.by_range.reduce((best, entry, k) => (
                (entry.rw ?? 0) > (data.by_range[best].rw ?? 0) ? k : best
            ), 0);
            if ((data.by_range[poorest].rw ?? 0) >= RELATIVE_MISFIT_FACTOR * data.rw && poorest !== worst) {
                const where = data.by_range[poorest][`${axis}_A`];
                checks.push({
                    level: 'note',
                    text: `Relative to its own signal the fit is poorest at ${axis} = ${where[0]}–${where[1]}${axis === 'r' ? ' Å' : ' Å⁻¹'} `
                        + `(Rw ${data.by_range[poorest].rw}, ${roundSig(data.by_range[poorest].rw / data.rw, 2)}× the overall), `
                        + 'though it holds little of the squared residual: weak signal, noise, or unreproduced sharp detail.'
                });
            }
            if (data.noise?.sigma) {
                const above = data.largest_residuals.filter((entry) => entry.above_noise);
                const { noise } = data;
                const against = noise.sigma === 'file'
                    ? 'the G(r) file\'s own σ'
                    : `σ(r) propagated from ${noise.from || 'F(Q)'}'s σ (RMCProfile's transform reproduces this G(r) to `
                        + `${noise.reproduction?.residual} of its maximum, tolerance ${noise.reproduction?.tolerance})`;
                const correlated = noise.correlation
                    ? `The r points are correlated (ρ = ${noise.correlation.neighbour} between neighbours, `
                        + `≈ ${noise.correlation.effective_points} independent of ${noise.correlation.points})`
                    : 'The r points of a Fourier-transformed G(r) are correlated';
                checks.push({
                    level: 'note',
                    text: `Against ${against}: ${above.length} of the ${data.largest_residuals.length} largest residual features `
                        + `${above.length === 1 ? 'is' : 'are'} above noise (≥ ${ABOVE_NOISE_SIGMA} σ`
                        + `${above.length ? `: ${axis} = ${above.map((entry) => entry[axis]).join(', ')}${unit}` : ''}); `
                        + `the residual's rms is ${noise.rms_residual_over_sigma} σ. ${correlated}: adjacent points beyond σ `
                        + 'are one feature, not several.'
                });
            } else if (data.noise) {
                checks.push({
                    level: 'note',
                    text: `No σ(r) for this G(r) (${data.noise.why}): whether a residual feature is above noise is not judged.`
                });
            }
            checks.unshift({ level: 'ok', text: `Rw ${data.rw} over ${axis} = ${data.range[0]}–${data.range[1]}${unit}; ${where}` });
            return { data, checks };
        }
    };
};

/** The tools over this run's props, in a stable order; absent data → no tool. */
export const buildRunTools = (run, { cards = [] } = {}) => {
    const context = buildRunContext(run || {});
    const catalog = Array.isArray(cards) ? cards : [];
    const cardIds = (test) => catalog.filter(test).map((card) => card.id);
    return [
        readRunDataTool(context, cardIds),
        convergenceTool(run?.rValueFile, cardIds),
        grPeaksTool(run, cardIds),
        fitTool(run?.plotFiles, cardIds)
    ].filter(Boolean);
};
