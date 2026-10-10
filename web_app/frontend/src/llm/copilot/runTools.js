// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

import { buildRunContext, seriesStats } from '../context/runContext';
import { extractPeaks } from '../context/pairCorrelations';
import { classifyConvergence } from '../watchdog/heuristics';

// Copilot tools over the run data this module already receives as props
// (plotFiles, rValueFile, structure, symmetry, runSettings, pcaSites) — so
// they need nothing from the host. The engine tools (sites, displacement
// directions, bond angles, symmetry at a tolerance) come from the host as the
// `tools` prop. Each tool is offered only when its data exists for this run.

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
    convergence: 'ln(chi^2) history of the last .log column'
};

// Least-squares slope per log row over the last `fraction` of the finite values.
const tailSlope = (values, fraction) => {
    const finite = values.filter(Number.isFinite);
    const count = Math.max(Math.min(finite.length, 3), Math.ceil(finite.length * fraction));
    const tail = finite.slice(-count);
    const n = tail.length;
    if (n < 2) return null;
    const meanX = (n - 1) / 2;
    const meanY = tail.reduce((sum, value) => sum + value, 0) / n;
    let covariance = 0;
    let variance = 0;
    tail.forEach((value, index) => {
        covariance += (index - meanX) * (value - meanY);
        variance += (index - meanX) ** 2;
    });
    return variance ? covariance / variance : 0;
};

const pairKey = (label) => {
    const match = /^([A-Z][a-z]?)\s*-\s*([A-Z][a-z]?)$/.exec((label || '').trim());
    return match ? [match[1], match[2]].sort().join('-') : null;
};

const readRunDataTool = (context) => {
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
        label: ({ section }) => `Read ${String(section || '').replace(/_/g, ' ')}`,
        run: async ({ section }) => ({ data: { section, [section]: context[section] }, checks: [] })
    };
};

const convergenceTool = (rValueFile) => {
    const history = rValueFile?.plotData?.series?.[0]?.y;
    if (!Array.isArray(history) || !history.length) return null;
    const column = rValueFile?.plotData?.chiColumn || null;
    return {
        name: 'assess_convergence',
        description: 'Classify the convergence of the ln(chi^2) history (the last .log column — one fit term, '
            + 'not the total) with the dashboard watchdog heuristics, and report the slope per log row over '
            + 'the last 10 %, 25 % and 50 % of the run. Status: improving | converged | stalled | diverging.',
        parameters: { type: 'object', properties: {} },
        label: () => 'Convergence check',
        page: 'dashboard',
        run: async () => {
            const stats = seriesStats(history);
            const status = classifyConvergence(history);
            const finite = history.filter(Number.isFinite);
            const data = {
                column,
                quantity: 'ln(chi^2), natural log; lower is better',
                status,
                n_steps: stats.nSteps,
                first: roundSig(stats.first),
                last: roundSig(stats.last),
                min: roundSig(stats.min),
                total_drop: finite.length >= 2 ? roundSig(finite[0] - finite[finite.length - 1]) : null,
                slope_per_step: {
                    last_10pct: roundSig(tailSlope(history, 0.1), 2),
                    last_25pct: roundSig(tailSlope(history, 0.25), 2),
                    last_50pct: roundSig(tailSlope(history, 0.5), 2)
                }
            };
            if (stats.nonFiniteSteps) data.non_finite_steps = stats.nonFiniteSteps;
            const checks = [{
                level: 'note',
                text: `Judged on ${column ? `the ${column} term` : 'the last .log column'} only, not the total chi^2.`
            }];
            if (status === 'unknown') checks.push({ level: 'fail', text: 'Fewer than two finite log rows — no trend to judge.' });
            else if (status === 'converged') checks.push({ level: 'ok', text: `Flat after a total drop of ${data.total_drop} in ln(chi^2).` });
            else if (status === 'improving') checks.push({ level: 'note', text: 'Still improving over the recent window.' });
            else checks.push({ level: 'warn', text: `Heuristic status: ${status}.` });
            if (stats.nonFiniteSteps) checks.push({ level: 'warn', text: `${stats.nonFiniteSteps} log rows are non-finite (NaN/Inf).` });
            if (stats.nSteps < 5) checks.push({ level: 'warn', text: `Only ${stats.nSteps} log rows.` });
            return { data, checks };
        }
    };
};

const grPeaksTool = (plotFiles, context) => {
    const partials = (plotFiles || []).find(
        (file) => file.plotKind === 'pdf_partials' && file.plotData?.series?.length
    );
    if (!partials) return null;
    const series = partials.plotData.series.filter((entry) => pairKey(entry.label));
    if (!series.length) return null;
    const labels = series.map((entry) => entry.label.trim());
    const constraints = context.run_settings?.min_distances_A;
    return {
        name: 'find_gr_peaks',
        description: 'Find the peaks of one partial g(r) (from the run\'s partial-PDF file) below r_max: '
            + 'position and approximate FWHM, a bond-length window around the first peak (use it for '
            + 'compute_bond_angles; a window over the first two peaks too when the second overlaps it), '
            + 'and the hard minimum-distance constraint of that pair when the run settings give one.',
        parameters: {
            type: 'object',
            properties: {
                pair: { type: 'string', enum: labels, description: 'Element pair, as labelled in the partials file.' },
                r_max: { type: 'number', minimum: 1, maximum: 20, description: 'Search up to this r (Å). Default 6.' },
                max_peaks: { type: 'integer', minimum: 1, maximum: 8, description: 'Default 3.' }
            },
            required: ['pair']
        },
        label: ({ pair }) => `g(r) peaks ${pair || ''}`.trim(),
        page: 'dashboard',
        run: async ({ pair, r_max: rMax = 6, max_peaks: maxPeaks = 3 }) => {
            const entry = series.find((item) => item.label.trim() === pair);
            const peaks = extractPeaks(entry.x, entry.y, { rMax, maxPeaks });
            const data = { pair, r_max_A: rMax, peaks_A: peaks.map((peak) => ({ r: roundSig(peak.r, 4), fwhm: roundSig(peak.fwhm, 2) })) };
            const checks = [];
            if (!peaks.length) {
                checks.push({ level: 'fail', text: `No peak above 15 % of the maximum below ${rMax} Å.` });
                return { data, checks };
            }
            const first = peaks[0];
            data.first_peak_window_A = [roundSig(Math.max(0, first.r - first.fwhm), 3), roundSig(first.r + first.fwhm, 3)];
            const key = pairKey(pair);
            const constraintEntry = constraints && !Array.isArray(constraints)
                ? Object.entries(constraints).find(([label]) => pairKey(label) === key)
                : null;
            if (constraintEntry && Number.isFinite(constraintEntry[1])) {
                const limit = constraintEntry[1];
                data.min_distance_constraint_A = limit;
                if (first.r - first.fwhm / 2 <= limit + 0.05) {
                    checks.push({
                        level: 'warn',
                        text: `The first peak's low side reaches the ${limit} Å closest-approach constraint — it may be constraint-limited.`
                    });
                }
            }
            // A second peak whose half-maximum reaches into the first-peak
            // window: the window cuts through it. Often both are one
            // coordination shell (a distorted polyhedron, e.g. the demo's
            // Ta–Se at 2.52 and 2.74 Å, where the first-peak window holds
            // only ~3.9 of the 6 Se), so offer a window over both.
            const second = peaks[1];
            if (second && second.r - second.fwhm / 2 < first.r + first.fwhm) {
                const reach = roundSig(second.r + second.fwhm / 2, 3);
                data.first_two_peaks_window_A = [data.first_peak_window_A[0], reach];
                checks.push({
                    level: 'note',
                    text: `The second peak (${roundSig(second.r, 4)} Å) overlaps the first-peak window, which cuts through it. `
                        + `If both are one coordination shell (a distorted polyhedron), use first_two_peaks_window_A `
                        + `(${data.first_peak_window_A[0]}–${reach} Å) and compare the coordination with the expected one.`
                });
            }
            if (!checks.length) checks.push({ level: 'ok', text: `First peak at ${roundSig(first.r, 4)} Å.` });
            return { data, checks };
        }
    };
};

/** The tools over this run's props, in a stable order; absent data → no tool. */
export const buildRunTools = (run) => {
    const context = buildRunContext(run || {});
    return [
        readRunDataTool(context),
        convergenceTool(run?.rValueFile),
        grPeaksTool(run?.plotFiles, context)
    ].filter(Boolean);
};
