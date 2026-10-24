// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// The Dashboard's G(r) error ribbon, the pure part: which σ an xPDF plot
// draws, the words of its caption and of its ? help (the tooltip with the
// transform, the reproduction residuals and their tolerance), and its Problems
// rows. components/GrRibbon.jsx draws the caption and the σ toggle.
//
// The σ comes from the plot's data as both runtimes give it: a series'
// `sigma` (the file's own σ column, browserData xpdfSigma / app._xpdf_sigma)
// and plotData.uncertainty (ftUncertaintyRun.js in the browser,
// /api/plot/data in Flask: rmc_toolkits/ft_uncertainty.py). With no σ in the
// file the ribbon is the propagated one; with the file's σ that one is drawn,
// and the propagated one is offered as a toggle.

import { rwpColumns } from './browserData';

export const RIBBON_LABELS = Object.freeze({ file: '±1σ (file)', propagated: '±1σ propagated' });

// The series the experiment's σ belongs to: the one carrying the file's σ,
// else the experimental column by the fit's roles (rwpColumns).
const experimentalSeries = (plotData) => {
    const series = plotData?.series || [];
    const withSigma = series.find((entry) => Array.isArray(entry.sigma) || ArrayBuffer.isView(entry.sigma));
    if (withSigma) return withSigma;
    const roles = rwpColumns(['x', ...series.map((entry) => entry.label)], series.length + 1);
    return roles ? series[roles[1] - 1] || null : null;
};

/**
 * What an xPDF plot can draw: { series, fileSigma, propagated, uncertainty },
 * or null for any other plot. `fileSigma` / `propagated` are σ arrays or null.
 */
export const ribbonSources = (file) => {
    const plotData = file?.plotData;
    if (file?.plotKind !== 'xpdf' || !plotData) return null;
    const series = experimentalSeries(plotData);
    const uncertainty = plotData.uncertainty || null;
    const propagated = uncertainty?.status === 'propagated' && Array.isArray(uncertainty.sigma)
        && series && uncertainty.sigma.length === series.x.length ? uncertainty.sigma : null;
    const fileSigma = series && (Array.isArray(series.sigma) || ArrayBuffer.isView(series.sigma)) ? series.sigma : null;
    return { series, fileSigma, propagated, uncertainty };
};

/** The σ the plot draws for a choice ('file' | 'propagated' | null): the file's by default. */
export const ribbonChoice = (sources, choice) => {
    if (!sources) return null;
    if (choice === 'propagated' && sources.propagated) return 'propagated';
    if (sources.fileSigma) return 'file';
    return sources.propagated ? 'propagated' : null;
};

/** InteractivePlot's `ribbon` prop for that choice, or null. */
export const ribbonProp = (sources, choice) => {
    const which = ribbonChoice(sources, choice);
    if (!which) return null;
    return {
        series: sources.series.label,
        sigma: which === 'file' ? sources.fileSigma : sources.propagated,
        label: RIBBON_LABELS[which]
    };
};

const SUPERSCRIPT = { '-': '⁻', 0: '⁰', 1: '¹', 2: '²', 3: '³', 4: '⁴', 5: '⁵', 6: '⁶', 7: '⁷', 8: '⁸', 9: '⁹' };

/** 3 significant digits in plain notation from 0.001 to 10⁴, else 1.5×10⁻¹⁴. */
export const sci = (value) => {
    if (!Number.isFinite(value)) return '—';
    if (value === 0) return '0';
    const abs = Math.abs(value);
    if (abs >= 1e-3 && abs < 1e4) return String(Number(value.toPrecision(3)));
    const [mantissa, exponent] = value.toExponential(1).split('e');
    return `${mantissa.replace(/\.0$/, '')}×10${[...String(Number(exponent))].map((ch) => SUPERSCRIPT[ch]).join('')}`;
};

const fixed = (value, digits = 2) => (Number.isFinite(value) ? value.toFixed(digits) : '—');
const capitalize = (text) => (text ? text[0].toUpperCase() + text.slice(1) : text);

/** The ? help: the transform, the reproduction tests with their tolerance, the correlation. */
export const ribbonDetails = (uncertainty) => {
    if (!uncertainty) return [];
    const lines = [];
    const { transform, reproduction, correlation, convolutionCheck, sources = {} } = uncertainty;
    if (transform) {
        lines.push(`RMCProfile's transform: G(r) = Σ ΔQ·Q·F(Q)·sin(Qr) over Q ${fixed(transform.qMin)}–${fixed(transform.qMax)} Å⁻¹ `
            + `(ΔQ ${sci(transform.dQ)}, ${transform.nQ} points), no 2/π, no ρ₀, no Lorch; r ${fixed(transform.rMin)}–${fixed(transform.rMax)} Å`
            + `${transform.convolve && transform.radius ? `; data box-convolved at R = ${fixed(transform.radius)} Å` : ''}.`);
    }
    const tolerance = reproduction?.tolerance;
    if (reproduction?.gr) {
        const { gr } = reproduction;
        lines.push(`Reproduction: ${sources.fq || 'the run’s F(Q)'} → ${sources.transform || 'the G(r) file'}: largest deviation `
            + `${sci(gr.residual)} of max |G| (tolerance ${sci(tolerance)}) ${gr.pass ? '✓' : '✗'}; fitted offset ${sci(gr.offset)}.`);
    }
    if (reproduction?.fq) {
        const { fq } = reproduction;
        lines.push(`Input: ${sources.input || 'the input file'} → ${sources.fq || 'the run’s F(Q)'} through the renormalisation `
            + `s = ${fixed(fq.scale, 4)}${fq.convolve ? ' and the box convolution' : ''}: largest deviation ${sci(fq.residual)} `
            + `(tolerance ${sci(tolerance)}) ${fq.pass ? '✓' : '✗'}.`);
    }
    if (correlation) {
        lines.push(`±1σ is the diagonal of the full covariance; the points are correlated: ρ = ${fixed(correlation.neighbour)} `
            + `between neighbours (Δr ${sci(correlation.gridStep)} Å), ${fixed(correlation.nyquist)} at π/Qmax = ${fixed(correlation.nyquistStep, 3)} Å; `
            + `≈ ${Math.round(correlation.effectivePoints)} independent points of ${correlation.points} `
            + `(Nyquist count ${Math.round(correlation.nyquistPoints)}).`);
    }
    if (convolutionCheck) {
        lines.push(`The box convolution, left out of the propagation, changes σ by ${sci(100 * convolutionCheck.maxChange)} % `
            + `at most (checked at ${convolutionCheck.points} points; limit ${sci(100 * convolutionCheck.limit)} %).`);
    }
    if (uncertainty.status === 'propagated') {
        lines.push('Statistical error of the measured F(Q) only: the reduction\'s systematic errors are not in its σ.');
    }
    return lines;
};

/**
 * The caption under an xPDF plot for a choice: { which, text, details } —
 * what the ribbon is (and that the points are correlated), or why there is
 * none; `details` are the ? help's lines.
 */
export const ribbonCaption = (sources, choice) => {
    const which = ribbonChoice(sources, choice);
    const uncertainty = sources?.uncertainty || null;
    let text;
    if (which === 'file') {
        text = 'Error ribbon: ±1σ from the file’s own σ column.';
    } else if (which === 'propagated') {
        const neighbour = uncertainty?.correlation?.neighbour;
        text = `Error ribbon: ±1σ propagated from ${uncertainty?.sources?.input || 'F(Q)'}’s σ through RMCProfile’s transform; `
            + `the r points are correlated${Number.isFinite(neighbour) ? ` (ρ = ${fixed(neighbour)} between neighbours)` : ''}.`;
    } else if (uncertainty?.reason) {
        text = `No error ribbon. ${capitalize(uncertainty.reason)}.`;
    } else {
        text = 'No error ribbon yet: reading the run’s F(Q) files…';
    }
    const details = ribbonDetails(uncertainty);
    if (which === 'file' && uncertainty?.reason) details.unshift(`No propagated ribbon. ${capitalize(uncertainty.reason)}.`);
    return { which: which || uncertainty?.status || 'pending', text, details };
};

/** Problems-section rows for the xPDF plots whose ribbon failed (a test, a check or a read). */
export const ribbonIssues = (files) => (files || [])
    .filter((file) => file.plotKind === 'xpdf' && file.plotData?.uncertainty?.status === 'failed')
    .map((file) => ({
        key: `ribbon:${file.path}`,
        message: `No error ribbon. ${capitalize(file.plotData.uncertainty.reason)}.`,
        severity: 'warning',
        source: file.name
    }));
