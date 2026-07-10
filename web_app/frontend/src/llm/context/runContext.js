// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// The seam where dashboard state becomes LLM input. buildRunContext distills
// the run object the dashboard already holds (parsed plot files, the combined
// R-value history, the folded structure, the detected symmetry) into a compact
// plain object whose JSON form fits comfortably in a small local model's
// context window. Everything here is a pure function of plain data — no React,
// no network — which is what makes it the most testable part of the module.

import { pairCorrelationsContext } from './pairCorrelations';

// Rough character budget for the serialized context: ~4.5k chars leaves even a
// small local model plenty of room for the system prompt, the question, and
// its answer, while fitting the symmetry/site and pair-correlation evidence.
export const CONTEXT_CHAR_BUDGET = 4500;

const HISTORY_POINTS = 48;
const MAX_SITES = 12;

const roundSig = (value, digits = 3) => (
    Number.isFinite(value) && value !== 0 ? Number(value.toPrecision(digits)) : value
);

// Uniform-stride downsample that always keeps the first and last points, with
// values rounded to 3 significant digits — raw float noise just wastes tokens.
export const downsampleSeries = (values, maxPoints = HISTORY_POINTS) => {
    if (!Array.isArray(values) || !values.length) return [];
    if (values.length <= maxPoints) return values.map((value) => roundSig(value));
    const picked = [];
    const step = (values.length - 1) / (maxPoints - 1);
    for (let i = 0; i < maxPoints; i += 1) {
        picked.push(roundSig(values[Math.round(i * step)]));
    }
    return picked;
};

// Least-squares slope (per step) over the tail of the series — the last 20% of
// points, at least 5 — so "still improving vs flat" reflects recent behavior.
export const recentSlope = (values) => {
    if (!Array.isArray(values) || values.length < 2) return 0;
    const count = Math.max(Math.min(values.length, 5), Math.ceil(values.length * 0.2));
    const tail = values.slice(-count);
    const n = tail.length;
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

// Summary statistics computed on the FULL series (before downsampling), so
// nothing important is lost to sampling.
export const seriesStats = (values) => {
    if (!Array.isArray(values) || !values.length) return null;
    return {
        nSteps: values.length,
        first: values[0],
        last: values[values.length - 1],
        min: Math.min(...values),
        max: Math.max(...values),
        recentSlopePerStep: recentSlope(values)
    };
};

const vectorLength = (vector) => Math.sqrt(vector.reduce((sum, value) => sum + value * value, 0));

const angleBetween = (a, b) => {
    const denominator = Math.max(vectorLength(a) * vectorLength(b), 1e-12);
    const cosine = a.reduce((sum, value, index) => sum + value * b[index], 0) / denominator;
    return Math.acos(Math.max(-1, Math.min(1, cosine))) * (180 / Math.PI);
};

// Cell lengths/angles are recomputed from the lattice vectors with the same
// math as ModelSummary.jsx — duplicated deliberately: importing from a host
// component would break this module's import boundary.
const structureContext = (structure) => {
    if (!structure?.latticeVectors || !structure?.supercell) return null;
    const cell = {
        cell_A: structure.latticeVectors.map((vector, index) => (
            roundSig(vectorLength(vector) / Math.max(structure.supercell[index], 1), 4)
        )),
        angles_deg: [
            angleBetween(structure.latticeVectors[1], structure.latticeVectors[2]),
            angleBetween(structure.latticeVectors[0], structure.latticeVectors[2]),
            angleBetween(structure.latticeVectors[0], structure.latticeVectors[1])
        ].map((angle) => roundSig(angle, 4)),
        supercell: structure.supercell,
        total_atoms: structure.totalAtoms
    };
    if (structure.elementCounts && Object.keys(structure.elementCounts).length) {
        cell.elements = structure.elementCounts;
    }
    return cell;
};

// The detected-symmetry block: space group, the tolerance ladder (distortion
// magnitude + symmetry character), and the Wyckoff orbits with per-orbit rms
// displacement aggregated from structure.basis[].dispA — sorted largest first
// so the sites participating most in local distortion lead, and budget
// trimming keeps them. Works with a minimal { spaceGroup } too (Flask mode
// sends no basis/orbits; those fields are simply absent).
const symmetryContext = (structure, symmetry) => {
    if (!symmetry?.spaceGroup) return null;
    const block = {
        space_group: `${symmetry.spaceGroup}`
            + `${symmetry.spaceGroupNumber ? ` (No. ${symmetry.spaceGroupNumber})` : ''}`
    };
    if (symmetry.pointGroup) block.point_group = symmetry.pointGroup;
    if (Number.isFinite(symmetry.nSpace)) block.n_ops = symmetry.nSpace;
    if (Number.isFinite(symmetry.toleranceA)) block.tolerance_A = roundSig(symmetry.toleranceA, 2);
    if (Number.isFinite(symmetry.maxResidual)) block.max_residual_A = roundSig(symmetry.maxResidual, 2);
    if (Array.isArray(symmetry.ladder) && symmetry.ladder.length > 1) {
        block.ladder = symmetry.ladder.map((brick) => ({
            sg: brick.spaceGroup,
            holds_A: [roundSig(brick.from, 2), roundSig(brick.to, 2)],
            n_ops: brick.nSpace
        }));
    }
    if (Array.isArray(symmetry.orbits) && symmetry.orbits.length) {
        const basis = structure?.basis;
        const sites = symmetry.orbits.map((orbit) => {
            const site = { element: orbit.element, multiplicity: orbit.size };
            if (orbit.wyckoff) site.wyckoff = `${orbit.size}${orbit.wyckoff}`;
            if (orbit.site) site.site_sym = orbit.site;
            if (Array.isArray(orbit.rep)) site.frac = orbit.rep.map((value) => roundSig(value, 3));
            const disps = (orbit.members || [])
                .map((index) => basis?.[index]?.dispA)
                .filter(Number.isFinite);
            if (disps.length) {
                site.mean_disp_A = roundSig(disps.reduce((sum, value) => sum + value, 0) / disps.length, 2);
                site.max_disp_A = roundSig(Math.max(...disps), 2);
            }
            return site;
        }).sort((a, b) => (b.mean_disp_A ?? -1) - (a.mean_disp_A ?? -1));
        block.sites = sites.slice(0, MAX_SITES);
        if (sites.length > MAX_SITES) block.sites_omitted = sites.length - MAX_SITES;
    }
    return block;
};

// Run-history counters from the .rmc6f header. The derived numbers are the
// useful ones: acceptance ratio, and accepted moves per atom — the standard
// gauge of whether the configuration has been sampled long enough.
const configurationOptimizationContext = (structure) => {
    const moves = structure?.moves;
    if (!moves) return null;
    const block = {};
    if (Number.isFinite(moves.generated)) block.moves_generated = moves.generated;
    if (Number.isFinite(moves.tried)) block.moves_tried = moves.tried;
    if (Number.isFinite(moves.accepted)) block.moves_accepted = moves.accepted;
    if (Number.isFinite(moves.accepted) && moves.tried > 0) {
        block.acceptance_ratio = roundSig(moves.accepted / moves.tried, 2);
    }
    if (Number.isFinite(moves.accepted) && structure.totalAtoms > 0) {
        block.accepted_moves_per_atom = roundSig(moves.accepted / structure.totalAtoms);
    }
    if (Number.isFinite(moves.accumulatedTimeS)) {
        block.accumulated_time_h = roundSig(moves.accumulatedTimeS / 3600);
    }
    return Object.keys(block).length ? block : null;
};

// Element-pair labels in RMCProfile's upper-triangle row-wise order (Ga-Ga,
// Ga-Nb, Ga-Se, Nb-Nb, …) — the order MINIMUM_DISTANCES values are given in,
// matching the partial-PDF column order.
const pairLabels = (atoms) => {
    const labels = [];
    for (let i = 0; i < atoms.length; i++) {
        for (let j = i; j < atoms.length; j++) labels.push(`${atoms[i]}-${atoms[j]}`);
    }
    return labels;
};

// The run-control settings parsed from <stem>.dat (browserData.parseRunSettings).
// Minimum distances are labeled per element pair when the counts line up, so
// the model can spot g(r) peaks pinned at a closest-approach constraint.
const runSettingsContext = (settings) => {
    if (!settings) return null;
    const block = {};
    for (const key of ['title', 'material', 'phase', 'temperature']) {
        if (settings[key]) block[key] = settings[key];
    }
    const atoms = settings.atoms || [];
    if (settings.minimumDistancesA?.length) {
        const labels = pairLabels(atoms);
        block.min_distances_A = labels.length === settings.minimumDistancesA.length
            ? Object.fromEntries(settings.minimumDistancesA.map((d, i) => [labels[i], d]))
            : settings.minimumDistancesA;
    }
    if (settings.maximumMovesA?.length) {
        block.max_move_A = atoms.length === settings.maximumMovesA.length
            ? Object.fromEntries(settings.maximumMovesA.map((d, i) => [atoms[i], d]))
            : settings.maximumMovesA;
    }
    if (settings.timeLimit) block.time_limit = settings.timeLimit;
    if (settings.savePeriod) block.save_period = settings.savePeriod;
    if (settings.weightOptimization) block.weight_optimization = true;
    if (settings.flags?.length) block.flags = settings.flags.slice(0, 8);
    if (settings.datasets?.length) {
        block.fitted_data = settings.datasets.map((dataset) => {
            const entry = { type: dataset.block.replace(/_DATA$/, '').toLowerCase() };
            if (dataset.file) entry.file = dataset.file;
            if (dataset.fit_type) entry.fit = dataset.fit_type;
            return entry;
        });
    }
    return Object.keys(block).length ? block : null;
};

const datasetContext = (plotFiles) => (
    (plotFiles || [])
        .filter((file) => file.plotKind && file.plotKind !== 'r_value')
        .map((file) => {
            const entry = { kind: file.plotKind };
            const title = file.plotData?.title;
            if (title && title !== file.plotKind) entry.title = title;
            const rwpValue = file.plotData?.metrics?.rwp;
            if (Number.isFinite(rwpValue)) entry.rwp = roundSig(rwpValue);
            return entry;
        })
);

const convergenceContext = (rValueFile, historyPoints) => {
    const history = rValueFile?.plotData?.series?.[0]?.y;
    const stats = seriesStats(history);
    if (!stats) return null;
    const convergence = {
        // The dashboard stores ln(chi^2), not raw chi^2 (browserData.js applies
        // Math.log when parsing the .log files) — say so, or the model will
        // misread the magnitudes.
        quantity: 'ln of chi^2 goodness metric (natural log; lower is better)',
        n_steps: stats.nSteps,
        first: roundSig(stats.first),
        last: roundSig(stats.last),
        min: roundSig(stats.min),
        max: roundSig(stats.max),
        recent_slope_per_step: roundSig(stats.recentSlopePerStep, 2),
        history: downsampleSeries(history, historyPoints)
    };
    const finalChi = rValueFile?.plotData?.metrics?.final_chi_r;
    if (Number.isFinite(finalChi)) convergence.final_chi_squared = roundSig(finalChi);
    return convergence;
};

// Distill dashboard props into the compact context object sent to the model.
// Every input is optional: missing pieces (no structure file, no logs) simply
// omit their section, and the prompt templates handle absence explicitly.
export const buildRunContext = ({
    runName,
    plotFiles = [],
    rValueFile = null,
    structure = null,
    symmetry = null,
    runSettings = null,
    liveData = false
} = {}) => {
    const context = {
        run: runName || 'unnamed run',
        live_mode: Boolean(liveData)
    };
    const structureInfo = structureContext(structure);
    if (structureInfo) context.structure = structureInfo;
    const optimization = configurationOptimizationContext(structure);
    if (optimization) context.configuration_optimization = optimization;
    const settingsInfo = runSettingsContext(runSettings);
    if (settingsInfo) context.run_settings = settingsInfo;
    const symmetryInfo = symmetryContext(structure, symmetry);
    if (symmetryInfo) context.symmetry = symmetryInfo;
    const pairs = pairCorrelationsContext(structure, plotFiles);
    if (pairs) context.pair_correlations = pairs;
    const datasets = datasetContext(plotFiles);
    if (datasets.length) context.datasets = datasets;
    const convergence = convergenceContext(rValueFile, HISTORY_POINTS);
    if (convergence) context.convergence = convergence;
    return context;
};

// Serialize with the character budget enforced, trimming the least essential
// detail first and never silently: extra g(r) peaks, then middle ladder rungs,
// then convergence-history points, then low-displacement sites, then the
// dataset list (each truncation recorded with an *_omitted count).
export const contextToJson = (context, budget = CONTEXT_CHAR_BUDGET) => {
    let current = context;
    let json = JSON.stringify(current, null, 1);
    const shrink = (patch) => {
        current = patch(current);
        json = JSON.stringify(current, null, 1);
    };

    if (json.length > budget && current.pair_correlations?.some((pair) => pair.gr_peaks_A?.length > 1)) {
        shrink((c) => ({
            ...c,
            pair_correlations: c.pair_correlations.map((pair) => (
                pair.gr_peaks_A?.length > 1 ? { ...pair, gr_peaks_A: pair.gr_peaks_A.slice(0, 1) } : pair
            ))
        }));
    }
    if (json.length > budget && current.symmetry?.ladder?.length > 2) {
        shrink((c) => ({
            ...c,
            symmetry: {
                ...c.symmetry,
                ladder: [c.symmetry.ladder[0], c.symmetry.ladder[c.symmetry.ladder.length - 1]],
                ladder_rungs_omitted: c.symmetry.ladder.length - 2
            }
        }));
    }
    for (const historyPoints of [24, 12]) {
        if (json.length <= budget) return json;
        if (current.convergence?.history?.length > historyPoints) {
            shrink((c) => ({
                ...c,
                convergence: {
                    ...c.convergence,
                    history: downsampleSeries(c.convergence.history, historyPoints)
                }
            }));
        }
    }
    if (json.length > budget && current.symmetry?.sites?.length > 8) {
        shrink((c) => ({
            ...c,
            symmetry: {
                ...c.symmetry,
                sites: c.symmetry.sites.slice(0, 8),
                sites_omitted: (c.symmetry.sites_omitted || 0) + c.symmetry.sites.length - 8
            }
        }));
    }
    if (json.length > budget && current.datasets?.length > 12) {
        shrink((c) => ({
            ...c,
            datasets: c.datasets.slice(0, 12),
            datasets_omitted: c.datasets.length - 12
        }));
    }
    return json;
};
