// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// The seam where dashboard state becomes LLM input. buildRunContext distills
// the run object the dashboard already holds (parsed plot files, the combined
// R-value history, the folded structure, the detected symmetry) into a compact
// plain object whose JSON form fits comfortably in a small local model's
// context window. Everything here is a pure function of plain data — no React,
// no network — which is what makes it the most testable part of the module.

import { pairCorrelationsContext } from './pairCorrelations';
import { labelledMinimumDistances } from './shells';
import { partialWeights } from './scatteringWeights';
// A cycle (heuristics imports seriesStats from here), safe because neither
// module calls the other while it is being evaluated.
import { classifyConvergence } from '../watchdog/heuristics';
import { describePointGroup } from './pointGroups';

// Rough character budget for the serialized context: ~6k chars (~1.8k tokens) leaves even a
// small local model plenty of room for the system prompt, the question, and
// its answer, while fitting the symmetry/site and pair-correlation evidence.
export const CONTEXT_CHAR_BUDGET = 6000;

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
// nothing important is lost to sampling. A log row whose chi^2 is non-finite
// (NaN / Inf / Fortran overflow — a blown-up run; JSON null from Flask) stays
// in the series: `first`/`last` are the raw end values (so a non-finite latest
// value is visible as such), min/max and the slope use the finite values only,
// and `nonFiniteSteps` counts the rest. One pass for min/max — a spread onto
// the argument stack throws past ~10^5 points.
export const seriesStats = (values) => {
    if (!Array.isArray(values) || !values.length) return null;
    const finite = values.filter(Number.isFinite);
    let min = Infinity;
    let max = -Infinity;
    finite.forEach((value) => {
        if (value < min) min = value;
        if (value > max) max = value;
    });
    return {
        nSteps: values.length,
        first: values[0],
        last: values[values.length - 1],
        min: finite.length ? min : NaN,
        max: finite.length ? max : NaN,
        recentSlopePerStep: recentSlope(finite),
        nonFiniteSteps: values.length - finite.length
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

/**
 * The tolerance (Å) below which the symmetry finder's verdicts are sampling
 * noise. Each average site is the mean of `count` box copies with 3-D rms
 * displacement dispA, so its position carries a 3-D standard error
 * dispA/√count; two symmetry-related means differ by √2 times that from noise
 * alone, and the finder takes the worst case over every site and operation —
 * so three times the largest such difference. On the demo (GaTa4Se8, 1000
 * copies, ~0.23 Å) this is 0.031 Å, the tolerance where F-43m first holds.
 * Null without a basis carrying counts and displacements.
 */
export const symmetryNoiseFloor = (basis) => {
    if (!Array.isArray(basis)) return null;
    let worst = 0;
    basis.forEach((site) => {
        if (Number.isFinite(site?.dispA) && site.count > 0) worst = Math.max(worst, site.dispA / Math.sqrt(site.count));
    });
    return worst > 0 ? 3 * Math.SQRT2 * worst : null;
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
    // A structure the finder did not analyse (too many sites, over the
    // operation budget) says why, instead of a bare 'not analysed'.
    if (symmetry.skipped && symmetry.reason) block.note = symmetry.reason;
    if (Array.isArray(symmetry.ladder) && symmetry.ladder.length > 1) {
        block.ladder = symmetry.ladder.map((brick) => ({
            sg: brick.spaceGroup,
            ...(brick.pointGroup ? { pg: brick.pointGroup } : {}),
            holds_A: [roundSig(brick.from, 2), roundSig(brick.to, 2)],
            n_ops: brick.nSpace
        }));
    }
    // What each point group implies (inversion, chirality, polarity), so an
    // answer never has to recall it (context/pointGroups.js).
    const groups = [symmetry.pointGroup, ...(symmetry.ladder || []).map((brick) => brick.pointGroup)];
    const facts = [...new Set(groups.filter(Boolean))].map(describePointGroup).filter(Boolean);
    if (facts.length) block.point_group_facts = facts;
    const floor = symmetryNoiseFloor(structure?.basis);
    if (floor) {
        block.noise_floor_A = roundSig(floor, 2);
        block.noise_floor_note = 'the site means carry sampling noise (rms displacement / √copies); a group change '
            + 'at a tolerance below noise_floor_A is within that noise, not evidence of a distortion';
    }
    if (Array.isArray(symmetry.orbits) && symmetry.orbits.length) {
        const basis = structure?.basis;
        const sites = symmetry.orbits.map((orbit) => {
            const site = { element: orbit.element, multiplicity: orbit.size };
            // `multiplicity` is the orbit's size in the GIVEN cell; the letter is
            // read in the standard cell the group is named in, so the label
            // pairs it with that cell's multiplicity (wyckoffMultiplicity) --
            // R3m's 3a for the Ga of a lacunar spinel kept in its F-cubic cell,
            // not 4a. The same rule as symmetryModel's orbitLabel(), which this
            // module may not import.
            if (orbit.wyckoff) site.wyckoff = `${orbit.wyckoffMultiplicity ?? orbit.size}${orbit.wyckoff}`;
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

// Per-site PCA of the RMC displacement clouds (the "PCA Ellipsoid" analysis):
// the anisotropic displacement parameters plus the non-Gaussianity that plain
// symmetry/mean-displacement stats can't express. `pcaSites` is the table the
// PCA Ellipsoid page computes (worker or Flask /api/pca/sites) — one entry per
// reference site with uIso, the three principal RMS amplitudes, anisotropy, and
// Mardia's multivariate excess kurtosis (b2 − 15)/5. Sites are ranked by
// |non-Gaussianity| (then uIso): a symmetric split site is NEGATIVE
// (platykurtic), so a signed ranking would list it last and trim it first.
const pcaContext = (pcaSites) => {
    const rows = pcaSites?.sites;
    if (!Array.isArray(rows) || !rows.length) return null;
    const sites = rows.map((site) => {
        const entry = {
            ref: site.referenceNumber,
            element: site.element,
            U_iso_A2: roundSig(site.uIso, 3),
            rms_axes_A: Array.isArray(site.rms) ? site.rms.map((value) => roundSig(value, 3)) : undefined,
            anisotropy: roundSig(site.anisotropy, 3)
        };
        if (Number.isFinite(site.nonGaussianity)) entry.non_gaussianity = roundSig(site.nonGaussianity, 3);
        if (site.degenerate) entry.degenerate = true;
        // No spread at all: anisotropy and non_gaussianity are undefined (null).
        if (site.zeroSpread) entry.zero_spread = true;
        // A mixed-occupancy site: `element` is the majority species only.
        if (site.mixed && site.elementCounts) {
            entry.mixed = true;
            entry.element_counts = site.elementCounts;
        }
        return entry;
    }).sort((a, b) => (
        Math.abs(b.non_gaussianity ?? 0) - Math.abs(a.non_gaussianity ?? 0)
        || (b.U_iso_A2 ?? -Infinity) - (a.U_iso_A2 ?? -Infinity)
    ));
    const block = {
        // Tell the model what the numbers mean, or it will misread the kurtosis.
        note: 'Per reference site, from PCA of RMC displacement clouds. U_iso_A2: isotropic '
            + 'ADP (Å²); rms_axes_A: principal RMS displacement amplitudes PC1≥PC2≥PC3 (Å); '
            + 'anisotropy = rms1/rms3; non_gaussianity: Mardia multivariate excess kurtosis '
            + '(b2-15)/5 — 0 = harmonic/Gaussian; >0 = peaked, heavy-tailed well or a minority '
            + 'off-centre component; <0 = flat-topped or bimodal, including a symmetric split site. '
            + 'Listed by |non_gaussianity|, largest first.',
        sites: sites.slice(0, MAX_SITES)
    };
    if (sites.length > MAX_SITES) block.sites_omitted = sites.length - MAX_SITES;
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

// The run-control settings parsed from <stem>.dat (browserData.parseRunSettings).
// Minimum distances are labeled per element pair when the counts line up (in
// RMCProfile's upper-triangle pair order, the labelling the shells' constraint
// reads: shells.js labelledMinimumDistances), so the model can spot g(r) peaks
// pinned at a closest-approach constraint.
const runSettingsContext = (settings) => {
    if (!settings) return null;
    const block = {};
    for (const key of ['title', 'material', 'phase', 'temperature']) {
        if (settings[key]) block[key] = settings[key];
    }
    const atoms = settings.atoms || [];
    const minDistances = labelledMinimumDistances(settings);
    if (minDistances) block.min_distances_A = minDistances;
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
            // One measurement fitted in both spaces is still one dataset:
            // its real-space G(r) is RMCProfile's transform of the same F(Q).
            if (dataset.real_space_fit && dataset.reciprocal_space_fit) {
                entry.fitted_in = 'reciprocal and real space (one measurement: the G(r) is its Fourier transform)';
            }
            return entry;
        });
        block.independent_datasets = settings.datasets.length;
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
            // <stem>_FT_XFQ1.csv is RMCProfile's Fourier transform of the
            // X-ray F(Q) it fits: the same measurement, not a second dataset.
            const name = file.name || file.path || '';
            const transformOf = /_FT_(XFQ\d+)\.csv$/i.exec(name);
            if (transformOf) entry.same_measurement_as = `${transformOf[1]} (its Fourier transform)`;
            return entry;
        })
);

// The Faber–Ziman weight of every partial in each radiation the run fits:
// what fraction of the measured signal a pair's correlations carry.
const partialWeightsContext = (structure, runSettings, plotFiles) => {
    const kinds = new Set((plotFiles || []).map((file) => file.plotKind));
    const blocks = (runSettings?.datasets || []).map((dataset) => dataset.block || '');
    const radiations = [];
    if (blocks.some((name) => /XRAY/.test(name)) || ['xpdf', 'xray_sq'].some((kind) => kinds.has(kind))) radiations.push('xray');
    if (blocks.some((name) => /NEUTRON/.test(name)) || ['npdf', 'neutron_sq', 'bragg'].some((kind) => kinds.has(kind))) radiations.push('neutron');
    const block = {};
    radiations.forEach((radiation) => {
        const weights = partialWeights(structure?.elementCounts, radiation);
        if (weights) block[radiation] = Object.fromEntries(weights.map(({ pair, weight }) => [pair, roundSig(weight, 2)]));
    });
    if (!Object.keys(block).length) return null;
    block.note = 'fraction of the measured signal each partial carries (Faber–Ziman; x-ray f≈Z at Q→0, '
        + 'neutron b_coh); a small weight means that pair is weakly constrained by the data';
    return block;
};

const convergenceContext = (rValueFile, historyPoints) => {
    const history = rValueFile?.plotData?.series?.[0]?.y;
    const stats = seriesStats(history);
    if (!stats) return null;
    // The series is the LAST column of the RMCProfile .log — the chi^2 of one
    // fit term, named by its header (e.g. X_ray_(R)1: the X-ray real-space
    // fit), not a total over datasets and constraints. Say which, or the model
    // will describe one term as the run's overall fit.
    const column = rValueFile?.plotData?.chiColumn || null;
    const convergence = {
        // The dashboard stores ln(chi^2), not raw chi^2 (browserData.js applies
        // Math.log when parsing the .log files) — say so, or the model will
        // misread the magnitudes.
        quantity: `ln of the chi^2 in the last .log column${column ? ` '${column}'` : ''} `
            + '(natural log; lower is better) — one fit term of the run, not a total',
        n_steps: stats.nSteps,
        first: roundSig(stats.first),
        last: roundSig(stats.last),
        min: roundSig(stats.min),
        max: roundSig(stats.max),
        // The dashboard watchdog's verdict on the same series (the source of
        // truth for "converged"), and what the recent window changed: a
        // model eyeballing 48 points reads a slow decline as a plateau.
        watchdog_status: classifyConvergence(history),
        recent_slope_per_step: roundSig(stats.recentSlopePerStep, 2),
        history: downsampleSeries(history, historyPoints)
    };
    const finite = history.filter(Number.isFinite);
    if (finite.length >= 10) {
        const tenth = Math.max(1, Math.round(finite.length * 0.1));
        convergence.change_over_last_10pct = roundSig(finite[finite.length - 1] - finite[finite.length - 1 - tenth], 2);
    }
    // Non-finite rows serialize as null; say how many and what they mean, or
    // the model reads a blown-up run's null tail as missing data.
    if (stats.nonFiniteSteps) {
        convergence.non_finite_steps = stats.nonFiniteSteps;
        convergence.non_finite_note = 'log rows whose chi^2 is NaN/Inf/overflow (null here) — '
            + 'the run produced non-finite values there; `last` is null when the latest row is one';
    }
    if (column) convergence.column = column;
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
    pcaSites = null,
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
    const pcaInfo = pcaContext(pcaSites);
    if (pcaInfo) context.pca_displacements = pcaInfo;
    const pairs = pairCorrelationsContext(structure, plotFiles, { runSettings });
    if (pairs) context.pair_correlations = pairs;
    const datasets = datasetContext(plotFiles);
    if (datasets.length) context.datasets = datasets;
    const weights = partialWeightsContext(structure, runSettings, plotFiles);
    if (weights) context.partial_weights = weights;
    const convergence = convergenceContext(rValueFile, HISTORY_POINTS);
    if (convergence) context.convergence = convergence;
    return context;
};

// Serialize with the character budget enforced, trimming the least essential
// detail first and never silently: extra g(r) peaks, then middle ladder rungs,
// then convergence-history points, then the least-anharmonic PCA sites (and, if
// still over, their explanatory note), then low-displacement symmetry sites,
// then the dataset list (each truncation recorded with an *_omitted count).
// JSON indented by one space, with each array of plain values, and then each
// object of plain values (and such arrays), on one line: the same text a model
// reads, at about half the length of one value per line. A bracket or brace
// followed by a line break never sits inside a JSON string (strings escape
// their line breaks), so only real arrays and objects are collapsed.
export const compactJson = (value) => JSON.stringify(value, null, 1)
    .replace(/\[\n\s*([^[\]{}]*?)\n\s*\]/g, (_, items) => `[${items.split(/,\n\s*/).join(', ')}]`)
    .replace(/\{\n\s*([^{}]*?)\n\s*\}/g, (_, items) => `{${items.split(/,\n\s*/).join(', ')}}`);

export const contextToJson = (context, budget = CONTEXT_CHAR_BUDGET) => {
    let current = context;
    let json = compactJson(current);
    const shrink = (patch) => {
        current = patch(current);
        json = compactJson(current);
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
    if (json.length > budget && current.pca_displacements?.sites?.length > 6) {
        shrink((c) => ({
            ...c,
            pca_displacements: {
                ...c.pca_displacements,
                sites: c.pca_displacements.sites.slice(0, 6),
                sites_omitted: (c.pca_displacements.sites_omitted || 0) + c.pca_displacements.sites.length - 6
            }
        }));
    }
    if (json.length > budget && current.pca_displacements?.note) {
        shrink((c) => {
            const { note, ...rest } = c.pca_displacements;
            return { ...c, pca_displacements: rest };
        });
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
