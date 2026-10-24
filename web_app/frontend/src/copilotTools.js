// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// The AI Copilot's analysis tools: the host side of the llm module's tool
// contract (src/llm/copilot/tools.js). Each tool runs one of the Workbench
// engines the pages use through the run's engine (`request`, the run module's
// `engine`: the Copilot's own PCA worker for a browser-loaded run or the Flask
// routes for a server folder, useCopilotRun.js; the workers' handlers
// in-thread in Node, openRunFolder.js), or the symmetry finder, and returns a
// compact summary for the model plus deterministic checks on the full result.
// The tools only compute: none changes a page, the shared tolerance, or a file.
// Each names the card it reads (`cards`), which breathes while it runs when
// that card is on screen (the llm module's card focus).
// The cluster and modulation engines decide what their statistics show (their
// `verdict`): analyze_clusters and analyze_modulations only shape the data and
// say each finding through a sentence table keyed by (finding, because).
// The summarizers are pure and exported for tests.

import { describeSymmetry, orbitLabel } from './symmetryModel';
import { siteLabel } from './siteLabel';
import { crystalOrientationRows, invert3, matVec3, transpose3, unitCellVectors } from './pcaCrystalFrame';
import { cnMatches, describePointGroup, lowerGroupTolerance, methodCheck, pointGroupFacts, shellCatalogue, windowCut } from './llm';
import { symmetryNoiseFloor } from './llm/context/runContext';
import { symmetryShifts } from './symmetryShift';
import { nearestLowIndex } from './viewer3d';
import { CLUSTER_BECAUSE, SPLIT_SITE_KURTOSIS, SPLIT_SITE_SE, isSplitSite } from './workers/clusters';
import { oneOf } from './llm/copilot/dataShape';

// Each tool's `dataShape` (llm/copilot/dataShape.js) is the shape of the data
// its summarizer builds — the paths the model quotes by reference — declared
// beside that summarizer and carried on the tool's spec.
const XYZ = ['number|null', 'number|null', 'number|null'];

const round = (value, digits = 3) => (
    Number.isFinite(value) ? (value === 0 ? 0 : Number(value.toPrecision(digits))) : null
);
const roundFixed = (value, decimals = 2) => (
    Number.isFinite(value) ? Number(value.toFixed(decimals)) + 0 : null
);

// Fewer displacement vectors than this make a sphere map mostly noise.
export const MIN_ORIENTATION_POINTS = 100;
// Calibrated one-sided sigma at which a peak or a whole-map departure counts.
export const SIGNIFICANCE_SIGMA = 3;

const parseWarningCheck = (payload) => (
    payload?.parseWarning ? [{ level: 'warn', text: `Structure file: ${payload.parseWarning}` }] : []
);

// --- split sites (list_sites) ----------------------------------------------------

// A site is split (two positions) by the app's one rule, isSplitSite in the
// cluster engine (workers/clusters.js, clusters.py's is_split_site): PC1
// resolved, its excess kurtosis ≤ SPLIT_SITE_KURTOSIS (−0.5) and at least
// SPLIT_SITE_SE (3) sampling errors (√(24/n)) below zero. analyze_clusters'
// split member sites are that engine's own verdict.
// The threshold as the texts the model reads write it (a true minus sign).
const SPLIT_RULE_KURTOSIS = String(SPLIT_SITE_KURTOSIS).replace('-', '−');

/**
 * A split site's static displacement along PC1 and its sampling error, from
 * its engine row (docs/algorithms/ai-assistant.md, Step 20):
 *
 *   δ² = λ₁ − (λ₂ + λ₃)/2 = rms1² − rms⊥²,
 *
 * rms⊥ the RMS of the two other axes: the thermal width along PC1 is taken as
 * the mean of the widths across it (exact for isotropic or uniaxial thermal
 * motion; with rms2 alone, the larger of two equal sample variances would
 * bias δ low). Each eigenvalue's sampling variance is (m4 − m2²)/n =
 * (κ + 2)·m2²/n, with the engine's population moments along that axis
 * (m2 = λ(n − 1)/n, κ its excess kurtosis): the error the engine's axis
 * resolution uses. The three add as independent (exact when the thermal
 * motion along PC1 is independent of the position taken and of the motion
 * across it), and the root by the delta method: SE(δ) = SE(δ²)/(2δ). For two
 * equally occupied positions ±δ this is the half distance between them; for
 * unequal ones, the RMS of the static displacement. Null without a spread
 * along PC1 beyond the other two axes.
 */
export const splitDisplacement = (site) => {
    const lambda = site?.eigenvalues;
    const n = site?.count;
    if (!Array.isArray(lambda) || lambda.length !== 3 || !lambda.every(Number.isFinite) || !(n > 1)) return null;
    const delta2 = lambda[0] - (lambda[1] + lambda[2]) / 2;
    if (!(delta2 > 0)) return null;
    // A collapsed axis has no kurtosis (κ null): its m2 is round-off, and the
    // Gaussian value adds nothing to the error.
    const variance = lambda.map((value, a) => {
        const kappa = site.excessKurtosis?.[a];
        const m2 = (value * (n - 1)) / n;
        return (Math.max((Number.isFinite(kappa) ? kappa : 0) + 2, 0) * m2 * m2) / n;
    });
    const delta = Math.sqrt(delta2);
    return { delta, se: Math.sqrt(variance[0] + (variance[1] + variance[2]) / 4) / (2 * delta) };
};

// --- list_sites ------------------------------------------------------------------

// The per-axis kurtosis only on resolved axes (an unresolved axis's direction
// is noise, and so is its κ), PC1's crystal direction only when PC1 is
// resolved, and delta_A / delta_se_A only on a split site.
const SITES_SHAPE = {
    note: 'string',
    'element?': 'string',
    n_sites: 'number',
    total_atoms: 'number|null',
    largest_U_iso_refs: ['number'],
    largest_abs_kurtosis_refs: ['number'],
    most_anisotropic_refs: ['number'],
    sites: [{
        ref: 'number',
        site: 'string',
        count: 'number',
        U_iso_A2: 'number|null',
        rms_axes_A: oneOf(XYZ, 'null'),
        anisotropy: 'number|null',
        kurtosis: 'number|null',
        'kurtosis_axes?': XYZ,
        'pc1_uvw?': ['number', 'number', 'number'],
        'delta_A?': 'number',
        'delta_se_A?': 'number',
        'zero_spread?': 'boolean'
    }]
};

// Split sites first, largest delta_A first; the rest by |kurtosis|, then
// U_iso — the order the run context uses — so trimming to the result budget
// drops the least remarkable sites first, and never a split site before
// another site.
const bySiteInterest = (a, b) => {
    if (a.delta_A != null || b.delta_A != null) {
        if (a.delta_A == null) return 1;
        if (b.delta_A == null) return -1;
        if (a.delta_A !== b.delta_A) return b.delta_A - a.delta_A;
    }
    return Math.abs(b.kurtosis ?? 0) - Math.abs(a.kurtosis ?? 0) || (b.U_iso_A2 ?? -Infinity) - (a.U_iso_A2 ?? -Infinity);
};

/**
 * The site table the model reads. `unitCell` (rows a, b, c in Å) gives PC1's
 * crystal direction; by default it is the payload's own (its supercell
 * lattice over the supercell).
 */
export const summarizeSites = (payload, { element = null, unitCell = undefined } = {}) => {
    const all = payload?.sites || [];
    const picked = element ? all.filter((site) => site.element === element || site.elementCounts?.[element]) : all;
    const cell = unitCell === undefined ? unitCellVectors(payload?.latticeVectors, payload?.supercell) : unitCell;
    const sites = picked.map((site) => {
        const row = {
            ref: site.referenceNumber,
            site: siteLabel(site),
            count: site.count,
            U_iso_A2: round(site.uIso),
            rms_axes_A: Array.isArray(site.rms) ? site.rms.map((value) => round(value)) : null,
            anisotropy: round(site.anisotropy),
            kurtosis: round(site.nonGaussianity)
        };
        const resolved = Array.isArray(site.axisResolved) ? site.axisResolved : [];
        const axes = [0, 1, 2].map((a) => (resolved[a] ? round(site.excessKurtosis?.[a]) : null));
        if (axes.some((value) => value != null)) row.kurtosis_axes = axes;
        // The sense whose closest crystal axis is acute (crystalOrientationRows),
        // as the PCA Ellipsoid table and analyze_displacement_directions quote it.
        const orientation = resolved[0] && Array.isArray(site.axes) && cell ? crystalOrientationRows(site.axes, cell) : null;
        if (orientation) row.pc1_uvw = orientation[0].crystalDirection.map((value) => roundFixed(value, 2));
        const split = isSplitSite(site) ? splitDisplacement(site) : null;
        if (split) {
            row.delta_A = round(split.delta);
            row.delta_se_A = round(split.se, 2);
        }
        if (site.zeroSpread) row.zero_spread = true;
        return row;
    }).sort(bySiteInterest);
    const ranked = (key) => [...sites]
        .filter((site) => Number.isFinite(site[key]))
        .sort((a, b) => Math.abs(b[key]) - Math.abs(a[key]))
        .slice(0, 3)
        .map((site) => site.ref);
    const splitRows = sites.filter((site) => site.delta_A != null);
    // The note defines only the fields this result holds.
    const note = [
        'kurtosis = Mardia multivariate excess kurtosis (0 Gaussian; >0 peaked/heavy-tailed; <0 flat or split)',
        sites.some((site) => site.kurtosis_axes)
            ? 'kurtosis_axes = excess kurtosis along PC1, PC2, PC3, null where that axis is not resolved (its direction is noise)'
            : null,
        'anisotropy = rms1/rms3',
        sites.some((site) => site.pc1_uvw) ? 'pc1_uvw = PC1 as a direct-lattice direction [u v w] (not a plane normal), where PC1 is resolved' : null,
        splitRows.length
            ? `delta_A = static displacement along PC1 of a split site (two positions: PC1 excess kurtosis ≤ ${SPLIT_RULE_KURTOSIS} and `
                + `${SPLIT_SITE_SE} sampling errors below 0), sqrt(rms1² − rms⊥²) with rms⊥² = (rms2² + rms3²)/2: ±delta_A for two `
                + 'equally occupied positions, their RMS otherwise; delta_se_A = its sampling error'
            : null,
        splitRows.length ? 'split sites first (largest delta_A first), then by |kurtosis|, largest first' : 'sites sorted by |kurtosis|, largest first'
    ].filter(Boolean).join('; ');
    const data = {
        note,
        ...(element ? { element } : {}),
        n_sites: sites.length,
        total_atoms: payload?.totalAtoms ?? null,
        largest_U_iso_refs: ranked('U_iso_A2'),
        largest_abs_kurtosis_refs: ranked('kurtosis'),
        most_anisotropic_refs: ranked('anisotropy'),
        sites
    };
    const checks = [...parseWarningCheck(payload)];
    if (!sites.length) {
        checks.push({
            level: 'fail',
            text: element ? `No site holds ${element}.` : 'The structure file has no reference sites.'
        });
    }
    if (payload?.reconstructed) {
        checks.push({ level: 'note', text: 'Sites were reconstructed by folding and clustering: the file has no site columns.' });
    }
    const flat = sites.filter((site) => site.zero_spread).length;
    if (flat) {
        checks.push({ level: 'warn', text: `${flat} site(s) have zero spread (an averaged configuration?): no axes or kurtosis.` });
    }
    if (splitRows.length) {
        const several = splitRows.length > 1;
        checks.push({
            level: 'note',
            text: `Split: site${several ? 's' : ''} ${splitRows.map((site) => `${site.ref} (${site.site})`).join(', ')} `
                + `${several ? 'have' : 'has'} two positions along PC1 (PC1 excess kurtosis ≤ ${SPLIT_RULE_KURTOSIS}, at least `
                + `${SPLIT_SITE_SE} sampling errors below 0); delta_A is the static displacement, along pc1_uvw.`
        });
    }
    if (sites.length && !checks.some((check) => check.level === 'warn' || check.level === 'fail')) {
        checks.push({ level: 'ok', text: `${sites.length} sites, ${payload?.totalAtoms ?? '?'} atoms.` });
    }
    // For the disorder check (the displacement-anisotropy interview,
    // 2026-10-24): which sites are split (the app's rule) or flat-topped (a
    // negative Mardia kurtosis that misses it), with their species.
    const facts = {
        sites: picked.map((site) => {
            const split = isSplitSite(site);
            return {
                ref: site.referenceNumber,
                element: site.element,
                elements: Object.keys(site.elementCounts || { [site.element]: 1 }),
                split,
                flat: !split && Number.isFinite(site.nonGaussianity) && site.nonGaussianity < 0
            };
        })
    };
    return { data, checks, facts };
};

// --- analyze_displacement_directions ---------------------------------------------

// A Cartesian direction as a direct-lattice [u v w] (fractional components
// M⁻ᵀ·d, largest magnitude 1), keeping its sign: the sphere map is never
// antipodally folded, so +d and −d are different directions here.
export const cartesianToUvw = (direction, unitCell) => {
    if (!Array.isArray(direction) || !unitCell) return null;
    const inverse = invert3(unitCell);
    if (!inverse) return null;
    const frac = matVec3(transpose3(inverse), direction);
    const largest = Math.max(...frac.map(Math.abs)) || 1;
    return frac.map((value) => roundFixed(value / largest, 2));
};

const anglesToCell = (direction, unitCell) => unitCell.map((edge) => {
    const length = Math.hypot(...edge) || 1;
    const norm = Math.hypot(...direction) || 1;
    const cosine = (direction[0] * edge[0] + direction[1] * edge[1] + direction[2] * edge[2]) / (length * norm);
    return roundFixed(Math.acos(Math.max(-1, Math.min(1, cosine))) * 180 / Math.PI, 1);
});

// The peak direction is the peak cell's centre: it names a direction only to
// within the cell's circumscribed radius (the engine's peakCellRadiusDeg, the
// largest angle from the centre to the cell's corners). Within
// LOW_INDEX_CONE_DEG the direction is resolved finely enough to tell the
// <100>, <110> and <111> families apart; a wider cell (every cell at ν ≤ 2,
// 19.3–37.4°) does not resolve it. Arm A+ applies the same bound
// (validation/armA.js, SWEEP_RULES.directionResolved).
export const directionResolved = (radiusDeg) => Number.isFinite(radiusDeg) && radiusDeg <= LOW_INDEX_CONE_DEG;

// --- a direction's tier (the displacement-anisotropy interview, 2026-10-24) -----

// The user's three tiers for a peak direction, on its unrounded σ and its
// enhancement: a finding at ≥ 3σ and ≥ 1.2× isotropic that holds at a
// neighbouring resolution; a hint at ≥ 3σ under 1.2×, or 2.5–3σ, or one
// that holds only at the engine's own resolution; nothing below 2.5σ.
export const DIRECTION_FINDING_ENHANCEMENT = 1.2;
export const DIRECTION_HINT_SIGMA = 2.5;

/**
 * The neighbouring resolution a finding is checked at: ν − 1, or ν + 1 at
 * ν = 1. The engine's resolution (recommended_frequency) is the finest tiling
 * whose cells keep at least 12 points each on average; ν − 1 keeps that
 * occupancy (coarser cells, more points each), ν + 1 would take it below the
 * engine's own floor.
 */
export const neighbourFrequency = (frequency) => (frequency > 1 ? frequency - 1 : frequency + 1);

// The angle (degrees) between two directions taken as axes (±u one axis).
const axisAngleDeg = (u, v) => {
    const norm = Math.hypot(...u) * Math.hypot(...v);
    if (!(norm > 0)) return null;
    const cosine = Math.abs(u[0] * v[0] + u[1] * v[1] + u[2] * v[2]) / norm;
    return Math.acos(Math.min(1, cosine)) * 180 / Math.PI;
};

/**
 * Whether a peak direction holds at the neighbouring resolution (`second`,
 * the engine's payload there): its own peak test passes (≥ 3σ, corrected over
 * its own cells) and its peak cell can hold the same axis as the first's —
 * the two directions, taken as axes, lie within the sum of the two cells'
 * radii (every direction assigned to a cell lies within its radius of its
 * centre). An axis, not a sense: the map is not antipodally folded, and a
 * site with two equal wells puts its peak on either well from one tiling to
 * the next (dev-jahn-teller's F: 148° apart at ν = 1 and 2); the +u/−u
 * imbalance has its own test. Returns { holds, apartDeg }.
 */
export const directionHolds = (first, second) => {
    const apartDeg = Array.isArray(first?.peakDirection) && Array.isArray(second?.peakDirection)
        ? axisAngleDeg(first.peakDirection, second.peakDirection) : null;
    const within = apartDeg != null && Number.isFinite(first.peakCellRadiusDeg) && Number.isFinite(second.peakCellRadiusDeg)
        && apartDeg <= first.peakCellRadiusDeg + second.peakCellRadiusDeg + 1e-9;
    return { holds: Boolean(within && second.peakSignificance >= SIGNIFICANCE_SIGMA), apartDeg };
};

/**
 * A peak direction's tier: { tier: 'finding' | 'hint' | 'none', because }.
 * `first` is the engine's payload at the call's resolution, `second` at the
 * neighbouring one (null when not asked), `auto` whether the call's
 * resolution is the engine's own: a finding is judged there.
 */
export const directionTier = (first, second = null, { auto = true } = {}) => {
    const sigma = first?.peakSignificance;
    if (!(sigma >= DIRECTION_HINT_SIGMA)) return { tier: 'none', because: 'below_hint_sigma' };
    if (sigma < SIGNIFICANCE_SIGMA) return { tier: 'hint', because: 'below_sigma' };
    if (!(first.peakEnhancement >= DIRECTION_FINDING_ENHANCEMENT)) return { tier: 'hint', because: 'small_enhancement' };
    if (!auto) return { tier: 'hint', because: 'not_auto' };
    if (!second) return { tier: 'hint', because: 'unchecked' };
    return directionHolds(first, second).holds ? { tier: 'finding', because: 'holds' } : { tier: 'hint', because: 'not_held' };
};

/** Whether a peak is a finding candidate: ≥ 3σ and ≥ 1.2× isotropic (what the neighbouring resolution then checks). */
export const directionCandidate = (payload) => payload?.peakSignificance >= SIGNIFICANCE_SIGMA
    && payload.peakEnhancement >= DIRECTION_FINDING_ENHANCEMENT;

// The peak direction as the checks write it: the crystal [u v w], else the
// Cartesian triple.
const peakDirectionText = (peak) => {
    if (Array.isArray(peak.direction_uvw)) return `[${peak.direction_uvw.join(' ')}]`;
    if (Array.isArray(peak.direction_cartesian)) return `(${peak.direction_cartesian.join(', ')}) (Cartesian)`;
    return 'the peak cell\'s centre';
};

// The engine's payload at the neighbouring resolution, as the result holds it.
const NEIGHBOUR_SHAPE = {
    frequency: 'number',
    significance_sigma: 'number|null',
    enhancement: 'number|null',
    radius_deg: 'number|null',
    apart_deg: 'number|null',
    holds: 'boolean'
};

// direction_uvw, angles_to_abc_deg and principal_axis_uvw need the unit cell;
// neighbour only when the neighbouring resolution was asked (a finding
// candidate at the engine's own resolution).
const ORIENTATION_SHAPE = {
    note: 'string',
    site: 'number|null',
    element: 'string',
    frequency: 'number',
    cells: 'number',
    used_points: 'number',
    rejected_points: 'number',
    enhancement_range: ['number|null', 'number|null'],
    peak: {
        enhancement: 'number|null',
        count: 'number',
        expected: 'number|null',
        significance_sigma: 'number|null',
        direction_cartesian: oneOf(XYZ, 'null'),
        'direction_uvw?': oneOf(XYZ, 'null'),
        'angles_to_abc_deg?': XYZ,
        radius_deg: 'number|null',
        tier: 'string'
    },
    'neighbour?': NEIGHBOUR_SHAPE,
    map_significance_sigma: 'number|null',
    antipodal_asymmetry: 'number|null',
    antipodal_asymmetry_z: 'number|null',
    antipodal_asymmetry_significant: 'boolean',
    orientation_tensor: {
        eigenvalues: ['number|null'],
        anisotropy: 'number|null',
        anisotropy_significance_sigma: 'number|null',
        'principal_axis_uvw?': XYZ
    },
    recommended_frequency: 'number'
};

/**
 * The directions result the model reads. `neighbour` is the engine's payload
 * at the neighbouring resolution (neighbourFrequency), asked for a finding
 * candidate at the engine's own resolution; `auto` whether the call ran at
 * the engine's own resolution. The peak's tier (directionTier) leads its
 * check, and the answer check reads it (facts.direction.tier).
 */
export const summarizeOrientation = (payload, unitCell, { pooledSites = 0, neighbour = null, auto = true } = {}) => {
    const direction = payload.peakDirection;
    const { tier, because } = directionTier(payload, neighbour, { auto });
    const peak = {
        enhancement: round(payload.peakEnhancement),
        count: payload.peakCount,
        expected: round(payload.peakExpected),
        significance_sigma: round(payload.peakSignificance),
        direction_cartesian: Array.isArray(direction) ? direction.map((value) => roundFixed(value, 3)) : null
    };
    if (unitCell && Array.isArray(direction)) {
        peak.direction_uvw = cartesianToUvw(direction, unitCell);
        peak.angles_to_abc_deg = anglesToCell(direction, unitCell);
    }
    // How finely the peak names its direction: the peak cell's radius.
    const radius = Number.isFinite(payload.peakCellRadiusDeg) ? roundFixed(payload.peakCellRadiusDeg, 1) : null;
    peak.radius_deg = radius;
    peak.tier = tier;
    const held = neighbour ? directionHolds(payload, neighbour) : null;
    const tensor = {
        eigenvalues: (payload.orientationEigenvalues || []).map((value) => round(value)),
        anisotropy: round(payload.orientationAnisotropy),
        anisotropy_significance_sigma: round(payload.orientationAnisotropySignificance)
    };
    const axisRows = unitCell && Array.isArray(payload.orientationAxes)
        ? crystalOrientationRows(payload.orientationAxes, unitCell)
        : null;
    if (axisRows) tensor.principal_axis_uvw = axisRows[0].crystalDirection.map((value) => roundFixed(value, 2));
    const data = {
        note: 'enhancement = density / isotropic (1 = isotropic); significances are calibrated one-sided sigma; '
            + 'the map is not antipodally folded, so +u and −u differ; [u v w] is a direct-lattice direction, not a plane normal; '
            + 'peak.radius_deg is the peak cell\'s radius: the peak names its direction only to within it; '
            + `peak.tier: finding (≥ ${SIGNIFICANCE_SIGMA}σ, ≥ ${DIRECTION_FINDING_ENHANCEMENT}×, held at the neighbouring resolution), hint, `
            + `or none (below ${DIRECTION_HINT_SIGMA}σ); the first check says why`
            + (neighbour ? '; neighbour: the peak at that resolution (held: ≥ 3σ, the same axis within the two cells\' radii)' : ''),
        site: payload.referenceNumber ?? null,
        element: payload.element ?? 'all',
        frequency: payload.frequency,
        cells: payload.cellCount,
        used_points: payload.usedPoints,
        rejected_points: payload.rejectedPoints,
        enhancement_range: [round(payload.vmin), round(payload.vmax)],
        peak,
        ...(neighbour ? {
            neighbour: {
                frequency: neighbour.frequency,
                significance_sigma: round(neighbour.peakSignificance),
                enhancement: round(neighbour.peakEnhancement),
                radius_deg: Number.isFinite(neighbour.peakCellRadiusDeg) ? roundFixed(neighbour.peakCellRadiusDeg, 1) : null,
                apart_deg: held.apartDeg == null ? null : roundFixed(held.apartDeg, 1),
                holds: held.holds
            }
        } : {}),
        map_significance_sigma: round(payload.mapSignificance),
        antipodal_asymmetry: round(payload.antipodalAsymmetry),
        antipodal_asymmetry_z: round(payload.antipodalAsymmetryZ),
        antipodal_asymmetry_significant: Boolean(payload.antipodalAsymmetrySignificant),
        orientation_tensor: tensor,
        recommended_frequency: payload.recommendedFrequency
    };
    const checks = [...parseWarningCheck(payload)];
    if (payload.usedPoints < MIN_ORIENTATION_POINTS) {
        checks.push({ level: 'warn', text: `Only ${payload.usedPoints} displacement vectors: the map is noisy.` });
    }
    const sigma = round(payload.peakSignificance);
    const strength = `the peak cell is ${round(payload.peakEnhancement)}× isotropic at ${sigma}σ (corrected over all cells)`;
    const where = `${peakDirectionText(peak)}${radius == null ? '' : ` ± ${radius}°, the cell radius at ν = ${payload.frequency}`}`;
    const atNeighbour = neighbour
        ? `at ν = ${neighbour.frequency} the peak is ${round(neighbour.peakSignificance)}σ, ${held.apartDeg == null ? 'its direction unknown' : `${roundFixed(held.apartDeg, 1)}° from it as an axis`} `
            + `(the two cells' radii ${roundFixed(payload.peakCellRadiusDeg + neighbour.peakCellRadiusDeg, 1)}°)`
        : '';
    if (tier === 'finding') {
        checks.push({
            level: 'ok',
            // A cell wider than LOW_INDEX_CONE_DEG cannot say which low-index
            // family the direction belongs to: say so, and name none.
            text: radius != null && !directionResolved(radius)
                ? `Preferred direction (a finding): ${strength}, toward ${where}; it holds: ${atNeighbour}. The direction is not resolved `
                    + `finer than ${radius}° (over ${LOW_INDEX_CONE_DEG}°), so it names no <100>, <110> or <111> family.`
                : `Preferred direction (a finding): ${where}; ${strength}; it holds: ${atNeighbour}.`
        });
    } else if (tier === 'hint') {
        const why = {
            below_sigma: `the peak is ${sigma}σ, between ${DIRECTION_HINT_SIGMA}σ and ${SIGNIFICANCE_SIGMA}σ`,
            small_enhancement: `${strength}, under the ${DIRECTION_FINDING_ENHANCEMENT}× a finding needs`,
            not_held: `${strength}, but it does not hold: ${atNeighbour}`,
            not_auto: `${strength}, but at ν = ${payload.frequency}, not the engine's resolution (ν = ${payload.recommendedFrequency}), where a finding is judged`,
            unchecked: `${strength}, not checked at a neighbouring resolution`
        }[because];
        checks.push({ level: 'note', text: `A hint of a direction, not a finding: ${why}; toward ${where}.` });
    } else {
        checks.push({
            level: 'note',
            text: `No preferred direction: the peak is ${sigma}σ, below ${DIRECTION_HINT_SIGMA}σ: nothing to report.`
        });
    }
    if (payload.mapSignificance == null) {
        checks.push({ level: 'note', text: 'The whole-map test is undefined here (too few points per cell).' });
    } else if (payload.mapSignificance >= SIGNIFICANCE_SIGMA) {
        // Any anisotropic Gaussian ellipsoid passes the map test (a
        // thermal-only O flattened across its bond: map 7.2σ, peak 2.3σ):
        // it shows anisotropy, and only a finding of the peak test a
        // preferred direction.
        checks.push({
            level: 'ok',
            text: `The displacements are anisotropic: the map departs from isotropy at ${round(payload.mapSignificance)}σ. `
                + (tier === 'finding'
                    ? (radius != null && !directionResolved(radius)
                        ? 'The peak test above places the preferred direction only to within its cell.'
                        : 'The peak test above names the preferred direction.')
                    : 'Thermal ellipsoids pass this test too: it is not a preferred direction, which needs a finding of the peak test.')
        });
    }
    if (payload.peakTieCount > 1) {
        checks.push({ level: 'note', text: `${payload.peakTieCount} cells tie for the peak: the peak direction is one of several equal ones.` });
    }
    if (payload.antipodalAsymmetrySignificant) {
        checks.push({
            level: 'note',
            text: `+u/−u imbalance at z = ${round(payload.antipodalAsymmetryZ)}: a skewness the ellipsoid cannot show.`
        });
    }
    // An element pooled over several sites can mix their local frames: each
    // vertex of a Ta4 tetrahedron sits on its own 3m axis, so a preference
    // every site shares in its own frame (along its own <111>) averages out.
    // The advice (which analysis to run instead) is method: arm C's method
    // 'off' sends the caveat alone (llm copilot/method.js).
    if (payload.referenceNumber == null && pooledSites > 1) {
        checks.push(methodCheck('note',
            `Pooled over ${pooledSites} sites. Where their local axes differ (the four vertices of a tetrahedral cluster), `
                + 'a direction each site prefers in its own frame cancels here, so an isotropic pooled map is not evidence '
                + 'against a per-site or cluster distortion.',
            { method: ' Test single sites (site_ref), or analyze_clusters for a cluster\'s shape.' }));
    }
    // For the answer check: no "preferred direction" unless one of the
    // question's analyses found one — a finding of the user's tiers
    // (2026-10-24); before, any peak at ≥ 3σ. The map test alone is
    // anisotropy (review L3, approved 2026-10-14).
    const subject = payload.referenceNumber != null ? `site ${payload.referenceNumber}` : (payload.element || 'all atoms');
    const facts = {
        direction: {
            subject,
            significant: payload.peakSignificance >= SIGNIFICANCE_SIGMA,
            anisotropic: payload.mapSignificance >= SIGNIFICANCE_SIGMA,
            peak_sigma: sigma,
            map_sigma: payload.mapSignificance == null ? null : round(payload.mapSignificance),
            tier,
            because,
            enhancement: round(payload.peakEnhancement),
            frequency: payload.frequency
        },
        // For the single-configuration check (the user's rule, 2026-10-24):
        // a distribution read from the run's one configuration.
        distribution: { kind: 'directions', subject: `the ${subject} displacement-direction distribution`, configurations: 1 }
    };
    return { data, checks, facts };
};

// --- compute_bond_angles ---------------------------------------------------------

// The strongest local maxima of a 5-bin running mean of the angle counts, each
// with the share of all angles within ±7.5° of it.
const angleModes = (centers, counts, total, maxModes = 3) => {
    if (!total) return [];
    const smooth = counts.map((_, index) => {
        let sum = 0;
        let n = 0;
        for (let k = -2; k <= 2; k += 1) {
            const j = index + k;
            if (j >= 0 && j < counts.length) { sum += counts[j]; n += 1; }
        }
        return sum / n;
    });
    const peaks = [];
    for (let i = 0; i < smooth.length; i += 1) {
        const left = i > 0 ? smooth[i - 1] : -Infinity;
        const right = i < smooth.length - 1 ? smooth[i + 1] : -Infinity;
        if (smooth[i] > 0 && smooth[i] >= left && smooth[i] > right) peaks.push(i);
    }
    return peaks
        .sort((a, b) => smooth[b] - smooth[a])
        .slice(0, maxModes)
        .map((index) => {
            const angle = centers[index];
            let inside = 0;
            centers.forEach((center, j) => { if (Math.abs(center - angle) <= 7.5) inside += counts[j]; });
            return { angle_deg: roundFixed(angle, 1), share: round(inside / total, 2) };
        })
        .sort((a, b) => a.angle_deg - b.angle_deg);
};

// The engine's counts re-binned into 18 bins of 10°, each bin's count shared
// by overlap: a 15° or 30° bin straddles two 10° bins, and putting it whole
// into the one its centre falls in left every third 10° bin empty (a comb).
export const tenDegreeBins = (centers, counts) => {
    const bins = new Array(18).fill(0);
    const width = centers.length > 1 ? centers[1] - centers[0] : 0;
    centers.forEach((center, index) => {
        const count = counts[index] || 0;
        if (!count) return;
        const lo = Math.max(0, center - width / 2);
        const hi = Math.min(180, center + width / 2);
        if (!(width > 0) || !(hi > lo)) {
            bins[Math.min(17, Math.max(0, Math.floor(center / 10)))] += count;
            return;
        }
        for (let j = Math.max(0, Math.floor(lo / 10)); j <= Math.min(17, Math.floor(hi / 10)); j += 1) {
            const overlap = Math.min(hi, 10 * (j + 1)) - Math.max(lo, 10 * j);
            if (overlap > 0) bins[j] += count * overlap / (hi - lo);
        }
    });
    return bins;
};

// The angle distribution as fractions in 18 bins of 10°.
const coarseHistogram = (centers, counts, total) => (
    tenDegreeBins(centers, counts).map((count) => (total ? round(count / total, 2) : 0))
);

// The same 18 bins divided by the share randomly oriented bonds would put
// there ((cos θ₁ − cos θ₂)/2): 1 everywhere for random directions. Raw counts
// deplete angles near 0° and 180° as sin θ, so an octahedron's trans angles
// hardly show in them; here they do.
export const sinCorrectedHistogram = (centers, counts, total) => {
    if (!total) return new Array(18).fill(0);
    const bins = tenDegreeBins(centers, counts);
    const rad = Math.PI / 180;
    return bins.map((count, j) => round((count / total) / ((Math.cos(10 * j * rad) - Math.cos(10 * (j + 1) * rad)) / 2), 2));
};

// For a central atom with about six or four neighbours (an A–B–A triplet in
// one window), the share of angles on each side of the split an ideal
// polyhedron has, against that ideal: an
// octahedron has 12 cis angles (90°) and 3 trans (180°) per centre, a
// tetrahedron six at 109.5°.
const POLYHEDRON_SPLIT_DEG = 135;
const polyhedronReference = (centers, counts, total, meanCoordination) => {
    if (!total || !Number.isFinite(meanCoordination)) return null;
    const band = (keep) => {
        let n = 0;
        let sum = 0;
        centers.forEach((center, index) => {
            if (keep(center)) { n += counts[index]; sum += counts[index] * center; }
        });
        return { share: round(n / total, 2), mean_deg: n ? roundFixed(sum / n, 1) : null };
    };
    if (cnMatches(meanCoordination, 6)) {
        return {
            ideal: 'octahedron: 80 % cis at 90°, 20 % trans at 180°',
            cis_below_135_deg: band((center) => center < POLYHEDRON_SPLIT_DEG),
            trans_above_135_deg: band((center) => center >= POLYHEDRON_SPLIT_DEG)
        };
    }
    if (cnMatches(meanCoordination, 4)) {
        return { ideal: 'tetrahedron: every angle at 109.5°', within_15_deg_of_109_5: band((center) => Math.abs(center - 109.47) <= 15) };
    }
    return null;
};

// The bonds' mean and spread; the spread from the engine's 40-bin length
// histogram over the window (bins ~0.01 Å for a 0.4 Å window), so the model
// quotes a measured width instead of inventing one.
const lengthSummary = (lengths) => {
    if (!lengths) return null;
    const summary = { mean_A: round(lengths.meanLength, 4), unique_bonds: lengths.uniqueBonds };
    const { binCenters, counts, count, meanLength } = lengths;
    if (count > 1 && Array.isArray(binCenters) && binCenters.length === counts?.length && Number.isFinite(meanLength)) {
        const variance = counts.reduce((sum, n, i) => sum + n * (binCenters[i] - meanLength) ** 2, 0) / count;
        summary.sd_A = round(Math.sqrt(variance), 2);
    }
    return summary;
};

const BOND_LENGTHS = oneOf({ mean_A: 'number|null', unique_bonds: 'number', 'sd_A?': 'number|null' }, 'null');
const ANGLE_BAND = { share: 'number|null', mean_deg: 'number|null' };
// A bond window against the run's first shell of its pair (the shell
// catalogue's window(): what it holds of the shell, its CN from the partial
// g(r) beside the shell's).
const SHELL_READING = { note: 'string', window_A: ['number', 'number'], cn: 'number|null', cn_in_window: 'number|null', holds: 'string' };
// A refused window keeps the counts and the shell readings that show the cut
// and drops the angles (mean_deg … histogram_10deg_sin_corrected), with
// `refused` saying why; the B–C window and bond read 'same as A–B' for one
// shared window; `shell_ab` only where the run reads a first shell for the
// A–B pair, `shell_bc` only for a B–C window of its own; `polyhedron` only for
// a whole shell of about six (octahedron) or four (tetrahedron).
const BOND_ANGLES_SHAPE = {
    triplet: 'string',
    window_ab_A: ['number', 'number'],
    window_bc_A: oneOf('string', ['number', 'number']),
    angle_count: 'number',
    central_atoms: 'number',
    'mean_deg?': 'number|null',
    'sd_deg?': 'number|null',
    'modes?': [{ angle_deg: 'number|null', share: 'number|null' }],
    'histogram_10deg?': ['number|null'],
    'histogram_10deg_sin_corrected?': ['number|null'],
    bond_ab: BOND_LENGTHS,
    bond_bc: oneOf('string', BOND_LENGTHS),
    coordination: {
        note: 'string', counts: ['number'], mean: 'number|null', modal: 'number|null', modal_share: 'number|null', unbonded_share: 'number|null'
    },
    'shell_ab?': SHELL_READING,
    'shell_bc?': SHELL_READING,
    'refused?': 'string',
    'polyhedron?': oneOf(
        { ideal: 'string', cis_below_135_deg: ANGLE_BAND, trans_above_135_deg: ANGLE_BAND },
        { ideal: 'string', within_15_deg_of_109_5: ANGLE_BAND }
    )
};

// A bond window against its first shell, as the data hold it (shell_ab, shell_bc).
const shellReadingOf = (reading, centre, partner) => ({
    note: `the first ${centre}–${partner} shell (find_gr_peaks); cn, cn_in_window: ${partner} per ${centre}`,
    window_A: [...reading.shell.window],
    cn: round(reading.shell.cn),
    cn_in_window: round(reading.cn),
    holds: reading.holds
});

/**
 * The bond-angle result the model reads, with deterministic checks. A window
 * that cuts the coordination shell is refused — its angles are not the
 * polyhedron's (the 2026-10-10 spot check read such a cut as "incomplete
 * octahedra") — unless the user chose it (`windowConfirmed`, the
 * coordination-shells skill's firm rule): the counts and shell readings that
 * show the cut are returned instead. The rule is the one the Workbench asks
 * on (llm/context/shells.js windowCut; approved by the user, 2026-10-19):
 * each bond window the tool reads (`shells`, the catalogue's window() for
 * each of its shellWindows: A–B, and B–C when it has a window of its own; null
 * where the run reads no shell) cuts its shell when it meets it and holds
 * fewer partner atoms than the whole shell by max(0.5, 15 % of the shell's
 * CN). Only where there is no shell to compare the A–B window with does the
 * engine's coordinationSummary decide (the counts fallback). The counts stay
 * in the data as information either way. For a whole shell of about six or
 * four neighbours the angles are also read against an ideal octahedron or
 * tetrahedron, and the sin θ-corrected histogram shows the trans angles raw
 * counts hide.
 */
export const summarizeBondAngles = (payload, { windowConfirmed = false, shells = [] } = {}) => {
    const total = payload.angleCount || 0;
    const coordination = Array.isArray(payload.coordination) ? payload.coordination : [];
    const apexTotal = coordination.reduce((sum, count) => sum + count, 0);
    // The engine's summary of the counts (triplets.py coordination_summary,
    // workers/triplets.js coordinationSummaryOf): the mean, the modal count
    // and its share, the share with no bond — information in the data, and
    // the counts fallback of the cut rule where the run reads no shell.
    const summary = payload.coordinationSummary;
    if (!summary || typeof summary !== 'object') {
        throw new Error('the bond-angle engine returned no coordinationSummary (triplets.py bond_angle_summary, workers/triplets.js bondAngleSummary)');
    }
    const meanCoordination = Number.isFinite(summary.mean) ? summary.mean : null;
    const [a, b, c] = payload.triplet || [];
    const data = {
        triplet: `${a}–${b}–${c} (${b} central)`,
        window_ab_A: payload.bond12,
        window_bc_A: payload.sharedEnds ? 'same as A–B' : payload.bond23,
        angle_count: total,
        central_atoms: payload.apexCount ?? apexTotal,
        mean_deg: roundFixed(payload.meanAngle, 1),
        sd_deg: roundFixed(payload.stdAngle, 1),
        modes: angleModes(payload.binCenters || [], payload.counts || [], total),
        histogram_10deg: coarseHistogram(payload.binCenters || [], payload.counts || [], total),
        histogram_10deg_sin_corrected: sinCorrectedHistogram(payload.binCenters || [], payload.counts || [], total),
        bond_ab: lengthSummary(payload.lengths12),
        bond_bc: payload.sharedEnds ? 'same as A–B' : lengthSummary(payload.lengths23),
        coordination: {
            note: 'counts[n] = central atoms with exactly n A–B bonds in the window',
            counts: coordination,
            mean: round(meanCoordination),
            modal: Number.isFinite(summary.modal) ? summary.modal : null,
            modal_share: round(summary.modalShare, 2),
            unbonded_share: round(summary.unbondedShare, 2)
        }
    };
    // Each bond window the tool reads (its shellWindows: A–B, and B–C when it
    // has a window of its own) by the one cut rule (shells.js windowCut, the
    // Workbench's): against the run's first shell of its pair, the A–B window
    // with the engine's counts as the fallback where there is none.
    const windows = [
        { key: 'shell_ab', partner: a, reading: shells[0] ?? null, counts: summary },
        ...(payload.sharedEnds ? [] : [{ key: 'shell_bc', partner: c, reading: shells[1] ?? null, counts: null }])
    ].map((item) => ({ ...item, verdict: windowCut(item.reading, item.counts) }));
    windows.forEach(({ key, partner, reading }) => {
        if (reading) data[key] = shellReadingOf(reading, b, partner);
    });
    const checks = [...parseWarningCheck(payload)];
    // Where the window comes from is method (llm copilot/method.js): arm C's
    // method 'off' sends the failure and the plain remedy.
    if (!total) {
        checks.push(methodCheck('fail',
            'No angles: no central atom has two neighbours in the window(s). Widen the window',
            { method: ' (find_gr_peaks gives the whole first shell, first_shell.window_A)' },
            ' or check the elements.'));
        return { data, checks };
    }
    // Why each cut window cuts: the CN shortfall against its shell (by path),
    // or, with no shell to compare, the counts (the fallback's words).
    const cutTexts = windows.filter(({ verdict }) => verdict.cut).map(({ key, partner, reading, verdict }) => {
        if (verdict.by === 'cn-margin') {
            return `The window cuts the first ${b}–${partner} shell: it holds ${round(reading.cn)} ${partner} per ${b} (${key}.cn_in_window) `
                + `of the shell's ${round(reading.shell.cn)} (${key}.cn, in ${key}.window_A), fewer by at least max(0.5, 15 % of the shell's CN).`;
        }
        if (verdict.by === 'unbonded') {
            return `${Math.round(summary.unbondedShare * 100)} % of the ${b} atoms have no ${a}–${b} bond in the window: it misses most of the shell.`;
        }
        const present = coordination.map((count, bonds) => (count > 0 ? bonds : null)).filter((bonds) => bonds != null);
        return `Only ${Math.round(summary.modalShare * 100)} % of the ${b} atoms have the most common count (${summary.modal} ${a}–${b} bonds); `
            + `the counts spread from ${present[0]} to ${present.at(-1)}. The window cuts through a coordination shell.`;
    });
    const cutText = cutTexts.length ? cutTexts.join(' ') : null;
    // For the coordination answer check (review L5): the A–B count in this
    // window, and whether it is the whole shell (its CN the shell's within the
    // CN margin; with no shell to compare, not cut by the counts) or a window
    // the user chose.
    const wholeShell = windows[0].verdict.whole;
    const facts = {
        coordination: {
            centre: b,
            partner: a,
            window: payload.bond12 ?? null,
            cn: round(meanCoordination),
            whole_shell: wholeShell,
            confirmed: windowConfirmed
        }
    };
    if (cutText && !windowConfirmed) {
        // Refused: keep what shows the cut, never the angles.
        ['mean_deg', 'sd_deg', 'modes', 'histogram_10deg', 'histogram_10deg_sin_corrected'].forEach((key) => { delete data[key]; });
        data.refused = 'the window cuts the coordination shell, so its angles are not the polyhedron\'s and are not returned';
        // Arm C's method 'off' still refuses; only the advice on which window goes.
        checks.push(methodCheck('fail',
            `Refused: ${cutText}`,
            { method: ' Use the whole shell (find_gr_peaks\' first_shell.window_A).' },
            ' If the user wants this window, ask them (ask_user) and call again with window_confirmed: true.'));
        return { data, checks, facts };
    }
    // For the single-configuration check (the user's rule, 2026-10-24): the
    // angle distribution and the coordination counts of the run's one
    // configuration (a refused window returns no distribution).
    if (total > 0) facts.distribution = { kind: 'bond_angles', subject: `the ${a}–${b}–${c} bond-angle distribution`, configurations: 1 };
    if (cutText) checks.push({ level: 'warn', text: `${cutText} Kept: the user chose this window; say that its angles are of a partial shell.` });
    if (!cutText) {
        checks.push({ level: 'ok', text: `${total} angles around ${apexTotal} ${b} atoms (mean coordination ${round(meanCoordination)}).` });
    }
    // Only a whole shell is read against an ideal polyhedron: a cut TaSe6
    // shell holds ~3.9 Se per Ta, which is not a tetrahedron. And only one
    // shell: an A–B–A triplet with one window. A mixed (Se–Ta–Te) or split
    // (two windows) triplet's angles are a subset of the polyhedron's.
    const oneShell = a === c && payload.sharedEnds === true;
    const reference = !wholeShell || !oneShell ? null : polyhedronReference(payload.binCenters || [], payload.counts || [], total, meanCoordination);
    if (reference) data.polyhedron = reference;
    if (reference?.trans_above_135_deg) {
        // A band with no angles has no mean: say "none", never "null°".
        const meanOf = (band) => (band.mean_deg == null ? 'none' : `${band.mean_deg}°`);
        // Which histogram to read the shape from is method (llm copilot/method.js):
        // arm C's method 'off' sends the counts against the ideal octahedron alone.
        checks.push(methodCheck('note',
            `Raw counts deplete angles near 180° (as sin θ): ${Math.round(reference.trans_above_135_deg.share * 100)} % of these angles `
                + `lie above 135° (mean ${meanOf(reference.trans_above_135_deg)}) and ${Math.round(reference.cis_below_135_deg.share * 100)} % `
                + `below (mean ${meanOf(reference.cis_below_135_deg)}), against an ideal octahedron's 20 % trans at 180° and 80 % cis `
                + 'at 90°.',
            { method: ' Read the shape from histogram_10deg_sin_corrected (1 = random directions), not from the raw modes.' }));
    }
    return { data, checks, facts };
};

/**
 * compute_bond_angles' bond windows (its `shellWindows`): the A–B window
 * around the central atom B, and the B–C window when the call gives one. The
 * Workbench asks about these before the tool runs (llm clarify.js), and the
 * tool refuses on the same windows by the same rule; `apply` puts a chosen
 * window back into the arguments, confirmed when it is the user's own.
 */
const bondShellWindows = (args) => [
    {
        centre: args.apex,
        partner: args.end1,
        window: [args.r12_min, args.r12_max],
        apply: ([lo, hi], { confirmed = false } = {}) => ({ ...args, r12_min: lo, r12_max: hi, ...(confirmed ? { window_confirmed: true } : {}) })
    },
    ...(args.r23_min != null && args.r23_max != null ? [{
        centre: args.apex,
        partner: args.end2,
        window: [args.r23_min, args.r23_max],
        apply: ([lo, hi], { confirmed = false } = {}) => ({ ...args, r23_min: lo, r23_max: hi, ...(confirmed ? { window_confirmed: true } : {}) })
    }] : [])
];

// --- analyze_clusters ------------------------------------------------------------

// analyze_clusters' own result budget (characters of the tool message).
export const CLUSTER_RESULT_BUDGET = 4500;

/**
 * analyze_clusters' method advice (llm copilot/method.js; arm C reads none):
 * where its window comes from, and which questions it answers — finite
 * clusters, never a corner-sharing polyhedron or another extended network,
 * which has no finite clusters for its null (the coordination-shells
 * interview, 2026-10-24, answer 7: a distortion index there). `polyhedra`:
 * the run offers polyhedral_distortion, which the advice then names.
 */
// (A constant: written inline before a ternary's colon, the tool name read as
// a result field to the skill lint's key scan.)
const POLYHEDRA_ADVICE = ': use polyhedral_distortion';
export const clusterMethod = ({ polyhedra = false } = {}) => [
    { after: 'in the average structure', text: ' (take the window from find_gr_peaks: first_shell.window_A of that element pair)' },
    'Use it for finite clusters (dimers, trimers, M4 tetrahedra, isolated polyhedra): their distortion, Jahn–Teller, '
        + 'dimerization or trimerization; single-site displacement directions cannot show these. Corner-sharing polyhedra '
        + `or another extended network have no finite clusters for its null${polyhedra ? POLYHEDRA_ADVICE : ''}.`
];

// One statistic against its null: change = obs/null − 1, when both are finite.
const NULL_ROW = { obs: 'number|null', null: 'number|null', z: 'number|null', 'change?': 'number' };
// No clusters (an extended network, or no two sites in the window): element,
// window_A, sites and sites_without_bond only; the rest comes together. The
// axis families are keyed <100>, <110>, <111>.
const CLUSTERS_SHAPE = {
    element: 'string',
    window_A: ['number', 'number'],
    sites: 'number',
    'sites_without_bond?': 'number',
    'clusters?': {
        members: 'number',
        edges: 'number',
        count: 'number',
        templates: 'number',
        'other_sizes?': [{ members: 'number', count: 'number' }]
    },
    'average_structure?': { edges_A: ['number|null'], spread_A: 'number|null', edge_se_A: 'number|null', spread_over_se: 'number|null' },
    'local?': {
        note: 'string',
        ladder_A: { obs: ['number|null'], null: ['number|null'], z: ['number|null'] },
        ladder_width_A: NULL_ROW,
        shape_rms_A: NULL_ROW,
        breathing_sd_A: NULL_ROW,
        mean_edge_A: NULL_ROW,
        'tetrahedral_modes?': { note: 'string', e_rms_A: NULL_ROW, t2_rms_A: NULL_ROW, t2_trigonality: NULL_ROW },
        'long_axis_within_15deg?': { '<*>': NULL_ROW },
        'short_axis_within_15deg?': { '<*>': NULL_ROW }
    },
    'outliers?': { bonds_outside_window: 'number|null', contacts_between_clusters_per_atom: 'number|null' },
    'split_member_sites?': ['number']
};

// How a favoured gyration axis is said, and the size the shape is read at.
const AXIS_VERB = { long: 'elongated', short: 'compressed' };
const SIZE_TEXT = '(at each cluster\'s own size)';

// The shape readings (one per cluster result: the edge spread against the
// shape null, with the pattern beside it). `width` is the data's rounded row.
const lessReading = ({ width, pattern }) => ({
    level: 'ok',
    text: `The clusters are ${Math.round(-width.change * 100)} % LESS distorted than independent atoms would make them ${SIZE_TEXT} `
        + `(edge spread ${width.obs} vs ${width.null} Å, z = ${width.z}): the atoms of a cluster move together`
        + (pattern
            ? `, yet ${pattern}: a correlated shape distortion whose amplitude the independent-atom null already holds (split or disordered sites).`
            : ', and there is no correlated shape distortion on top.')
});
const patternReading = (context) => (context.less ? lessReading(context) : {
    level: 'ok',
    text: `The edge spread (${context.width.obs} Å) matches the shuffle null ${SIZE_TEXT} (${context.width.null} Å; change `
        + `${Math.round((context.width.change ?? 0) * 100)} %, z = ${context.width.z}), but the pattern does not: ${context.pattern}. A correlated shape `
        + 'distortion whose amplitude the sites\' own displacements already hold.'
});
const noDistortion = ({ width }) => ({
    level: 'ok',
    text: `No cluster shape distortion beyond independent atoms ${SIZE_TEXT}: the edge spread (${width.obs} Å) matches the shuffle null `
        + `(${width.null} Å; change ${Math.round((width.change ?? 0) * 100)} %, z = ${width.z}).`
});
// A cluster with one edge (a dimer) has no edge spread: its shape reading, whatever `distorted` says.
const singleEdgeReading = ({ kind, pattern }) => ({
    level: 'ok',
    text: `${kind.members === 2 ? 'Dimers have' : 'Each cluster has'} one edge: the shape is the bond length, read from the size `
        + `row (breathing_sd_A, against independent atoms)${pattern ? `; ${pattern}` : ''}.`
});

/**
 * How analyze_clusters says each finding of the cluster engine's verdict
 * (workers/clusters.js, clusters.py): (finding, because) → a check, or null
 * where the finding adds no check of its own (rigid and lessDistorted are said
 * by the shape reading and the breathing note, favouredAxes by its families,
 * the negatives by nothing). Every `because` the engine can give has an entry,
 * so a new one cannot go unsaid (copilotTools.test.js checks). The texts are
 * the tool's from before the rules moved into the engines (2026-10-19).
 */
export const CLUSTER_SENTENCES = {
    // Which window, and where to take it from, is method (llm
    // copilot/method.js): arm C's method 'off' sends the fact and the plain remedy.
    'window.extendedNetwork': {
        extended_network: ({ element, lo, hi }) => methodCheck('fail',
            `The ${lo}–${hi} Å window links the ${element} atoms into an extended network (chains or a 3-D net), `
                + 'not finite clusters: narrow it',
            { method: ` to the intra-cluster bond (find_gr_peaks: first_shell.window_A of ${element}-${element})` },
            '.'),
        no_clusters: ({ element, lo, hi }) => methodCheck('fail',
            `No two ${element} sites of the average structure lie ${lo}–${hi} Å apart: no clusters. Widen the window`,
            { method: ` (find_gr_peaks gives the ${element}-${element} first shell)` },
            '.'),
        finite_clusters: null
    },
    'window.mixedSizes': {
        above_bar: ({ payload }) => ({
            level: 'warn',
            text: `The window gives clusters of different sizes (${payload.kinds.map((k) => `${k.clusters}×${k.members}`).join(', ')}): `
                + 'it may cut through or merge clusters. The statistics are for the most common size only.'
        }),
        within_bar: null,
        extended_network: null,
        no_clusters: null
    },
    staticDistortion: {
        above_sampling_error: ({ ideal }) => ({
            level: 'note',
            text: `The average cluster is itself distorted: its edges span ${ideal.edges_A[0]}–${ideal.edges_A.at(-1)} Å, `
                + `${ideal.spread_over_se}× their sampling error (${ideal.edge_se_A} Å) — a static distortion of the average structure.`
        }),
        within_sampling_error: ({ ideal }) => ({
            level: 'note',
            text: `The average cluster is regular within its sampling error: edges ${ideal.edges_A[0]}–${ideal.edges_A.at(-1)} Å `
                + `(± ${ideal.edge_se_A} Å each).`
        }),
        extended_network: null,
        no_clusters: null
    },
    distorted: {
        edge_spread: ({ width }) => ({
            level: 'ok',
            text: `The clusters are distorted beyond independent atoms ${SIZE_TEXT}: the spread of their edges (longest − shortest, `
                + `${width.obs} Å) is ${Math.round(width.change * 100)} % above the shuffle null (${width.null} Å, z = ${width.z}) — `
                + 'a correlated distortion the data put in.'
        }),
        axis_family: patternReading,
        trigonality: patternReading,
        less_distorted: ({ width }) => lessReading({ width, pattern: null }),
        within_null: noDistortion,
        undefined_statistic: noDistortion,
        extended_network: null,
        no_clusters: null
    },
    lessDistorted: Object.fromEntries(CLUSTER_BECAUSE.lessDistorted.map((because) => [because, null])),
    breathing: {
        above_null: ({ breathing }) => ({
            level: 'note',
            text: `The clusters breathe: their size (mean edge) fluctuates ${Math.round(breathing.change * 100)} % more than for independent `
                + `atoms (breathing sd ${breathing.obs} vs ${breathing.null} Å, z = ${breathing.z}) — a collective size fluctuation, `
                + 'not a shape distortion.'
        }),
        below_null: ({ breathing, split }) => ({
            level: 'note',
            text: `The cluster size (mean edge) fluctuates ${Math.round(-breathing.change * 100)} % less than for independent atoms `
                + `(breathing sd ${breathing.obs} vs ${breathing.null} Å, z = ${breathing.z})`
                + (split
                    ? ': the split sites\' positions combine in a pattern that keeps the size (as a Jahn–Teller distortion does), '
                        + 'not stiff bonds.'
                    : ': the bonds are stiffer than independent atoms make them — the atoms move together.')
        }),
        within_null: null,
        undefined_statistic: null,
        extended_network: null,
        no_clusters: null
    },
    rigid: Object.fromEntries(CLUSTER_BECAUSE.rigid.map((because) => [because, null])),
    trigonalLean: {
        leans_trigonal: ({ trig }) => ({
            level: 'note',
            text: `The edge pattern leans trigonal: three edges at one vertex differ from the other three (C3v-like; trigonality ${trig.obs} vs ${trig.null}, z = ${trig.z}).`
        }),
        leans_one_edge: ({ trig }) => ({
            level: 'note',
            text: `The edge pattern leans to one edge against its opposite (trigonality ${trig.obs} vs ${trig.null}, z = ${trig.z}).`
        }),
        within_null: null,
        not_tetrahedra: null,
        undefined_statistic: null,
        extended_network: null,
        no_clusters: null
    },
    favouredAxes: Object.fromEntries(CLUSTER_BECAUSE.favouredAxes.map((because) => [because, null])),
    'favouredAxes.families': {
        enhanced: ({ family, entry, verb }) => ({
            level: 'note',
            text: `More clusters are ${verb} along ${family} than by chance: ${Math.round(entry.obs * 100)} % `
                + `vs ${Math.round(entry.null * 100)} % in the null (z = ${entry.z}).`
        }),
        within_null: null,
        undefined_statistic: null
    },
    splitMemberSites: {
        // What to read and run next is method: arm C's method 'off' sends the caveat alone.
        half_or_more_split: ({ verdict }) => methodCheck('note',
            `Split sites: ${verdict.splitMemberSites.sites.length} of the ${verdict.splitMemberSites.memberSites.length} member sites `
                + `(${verdict.splitMemberSites.sites.join(', ')}) `
                + 'have two positions (PC1 excess kurtosis ≤ −0.5). The null keeps each site\'s two positions, so it already holds the '
                + 'distortion\'s amplitude: the edge spread cannot show it, and "no distortion" or "less distorted" does not make the '
                + 'clusters regular.',
            { method: ' Read the pattern (axis families, trigonality), the partial g(r) shell (find_gr_peaks) and the bond angles.' }),
        fewer_than_half_split: null,
        none_split: null,
        extended_network: null,
        no_clusters: null
    },
    'window.bondsOutsideWindow': {
        above_bar: ({ data, lo, hi }) => ({
            level: 'warn',
            text: `${Math.round(data.outliers.bonds_outside_window * 100)} % of the cluster bonds lie outside ${lo}–${hi} Å in the configuration: `
                + 'the bond lengths spread wider than the window, or have heavy tails (atoms straying from their cluster).'
        }),
        within_bar: null,
        extended_network: null,
        no_clusters: null
    },
    'window.contactsBetweenClusters': {
        above_bar: ({ data, element }) => ({
            level: 'warn',
            text: `${data.outliers.contacts_between_clusters_per_atom} ${element}–${element} contacts per atom inside the window `
                + 'are not cluster bonds: atoms crossing the gap between clusters.'
        }),
        within_bar: null,
        extended_network: null,
        no_clusters: null
    }
};

// One finding said into `checks` through a sentence table: an unknown
// `because` is a broken contract with the engine, never silent.
const sayInto = (checks, table, tool) => (name, finding, context) => {
    const entries = table[name];
    if (!entries || !Object.prototype.hasOwnProperty.call(entries, finding?.because)) {
        throw new Error(`${tool}: the engine's ${name} finding gave '${finding?.because}', which has no sentence`);
    }
    const check = entries[finding.because]?.(context);
    if (check) checks.push(check);
};

/**
 * The cluster result the model reads: the engine's statistics shaped into
 * rows (the shape rows against the engine's shape null — independent atoms at
 * each cluster's own breathing — and the size rows against independent atoms;
 * the review's L2, approved 2026-10-14) and its verdict said through
 * CLUSTER_SENTENCES. The engine decides; this only renders.
 */
export const summarizeClusters = (payload) => {
    const [lo, hi] = payload.window || [];
    const element = payload.element;
    const { verdict } = payload;
    const data = {
        element,
        window_A: [lo, hi],
        sites: payload.sites?.length ?? 0
    };
    if (payload.isolatedSites) data.sites_without_bond = payload.isolatedSites;
    const checks = [...parseWarningCheck(payload)];
    const say = sayInto(checks, CLUSTER_SENTENCES, 'analyze_clusters');
    if (verdict.window.extendedNetwork.because !== 'finite_clusters') {
        say('window.extendedNetwork', verdict.window.extendedNetwork, { element, lo, hi });
        return { data, checks };
    }
    const { observed: obs } = payload;
    const shape = { mean: payload.shapeNull.mean, z: payload.shapeNull.z };
    const size = { mean: payload.nullMean, z: payload.z };
    const row = (value, nullValue, zValue, digits = 3) => {
        const entry = { obs: round(value, digits), null: round(nullValue, digits), z: round(zValue, 2) };
        if (Number.isFinite(value) && Number.isFinite(nullValue) && nullValue !== 0) entry.change = round(value / nullValue - 1, 2);
        return entry;
    };
    const kind = payload.kind;
    data.clusters = {
        members: kind.members,
        edges: kind.edges,
        count: kind.clusters,
        templates: kind.templates
    };
    // The other kinds, never the one analysed (a tie once listed it here).
    const others = payload.kinds.filter((k) => k.members !== kind.members || k.edges !== kind.edges);
    if (others.length) data.clusters.other_sizes = others.map((k) => ({ members: k.members, count: k.clusters }));
    data.average_structure = {
        edges_A: payload.ideal.ladder.map((value) => round(value, 4)),
        spread_A: round(payload.ideal.spread, 2),
        edge_se_A: round(payload.ideal.edgeSe, 2),
        spread_over_se: round(payload.ideal.spreadOverSe, 2)
    };
    data.local = {
        note: 'obs = the configuration. Shape rows (ladder, edge spread, shape rms, modes, axis families): null = each '
            + 'cluster keeps its own breathing (size) and each atom its own site\'s other displacements, taken from a '
            + 'random other cluster. Size rows (breathing_sd_A, mean_edge_A): null = independent atoms (each atom\'s whole '
            + 'displacement from a random other cluster). z = (obs − null) / null sd over the shuffles; change = obs/null − 1; '
            + 'ladder_A[k] = the mean over clusters of each cluster\'s k-th shortest edge (not the range of single bonds)',
        ladder_A: {
            obs: obs.ladder.map((value) => round(value, 4)),
            null: shape.mean.ladder.map((value) => round(value, 4)),
            z: shape.z.ladder.map((value) => round(value, 2))
        },
        ladder_width_A: row(obs.ladderWidth, shape.mean.ladderWidth, shape.z.ladderWidth),
        shape_rms_A: row(obs.shapeRms, shape.mean.shapeRms, shape.z.shapeRms),
        breathing_sd_A: row(obs.breathingSd, size.mean.breathingSd, size.z.breathingSd),
        mean_edge_A: row(obs.meanEdge, size.mean.meanEdge, size.z.meanEdge, 4)
    };
    if (Number.isFinite(obs.eRms)) {
        data.local.tetrahedral_modes = {
            note: 'E: opposite-edge-pair sums (tetragonal-type); T2: opposite-pair differences; t2_trigonality = '
                + 'mean Σd⁴/(Σd²)², 1/3 when three edges at one vertex differ alike (trigonal, C3v), 1 when one edge differs from its opposite',
            e_rms_A: row(obs.eRms, shape.mean.eRms, shape.z.eRms),
            t2_rms_A: row(obs.t2Rms, shape.mean.t2Rms, shape.z.t2Rms),
            t2_trigonality: row(obs.t2Trigonality, shape.mean.t2Trigonality, shape.z.t2Trigonality)
        };
    }
    const familyRows = (name) => (obs[name] ? Object.fromEntries(Object.keys(obs[name]).map((family) => [
        `<${family}>`, row(obs[name][family], shape.mean[name][family], shape.z[name][family], 3)
    ])) : undefined);
    if (obs.longAxis) data.local.long_axis_within_15deg = familyRows('longAxis');
    if (obs.shortAxis) data.local.short_axis_within_15deg = familyRows('shortAxis');
    data.outliers = {
        bonds_outside_window: round(payload.outliers.brokenEdgeFraction, 2),
        contacts_between_clusters_per_atom: round(payload.outliers.extraContactsPerAtom, 2)
    };
    // The analysed kind's member sites the engine found split (all of them; the
    // split finding needs half).
    if (verdict.splitMemberSites.sites.length) data.split_member_sites = [...verdict.splitMemberSites.sites];

    // The pattern beside the shape reading: the first favoured axis family, or the trigonal lean.
    const favoured = verdict.favouredAxes.families.filter((family) => family.answer === true);
    const pattern = favoured.length
        ? `more clusters than by chance are ${AXIS_VERB[favoured[0].axis]} along <${favoured[0].family}> (below)`
        : (verdict.trigonalLean.answer === true ? 'the edge pattern leans to one form (below)' : null);
    const context = {
        element, lo, hi, payload, data, verdict, kind, pattern,
        width: data.local.ladder_width_A,
        breathing: data.local.breathing_sd_A,
        trig: data.local.tetrahedral_modes?.t2_trigonality,
        ideal: data.average_structure,
        less: verdict.lessDistorted.answer === true,
        split: verdict.splitMemberSites.answer === true
    };
    say('window.mixedSizes', verdict.window.mixedSizes, context);
    say('staticDistortion', verdict.staticDistortion, context);
    if (kind.edges === 1) checks.push(singleEdgeReading(context));
    else say('distorted', verdict.distorted, context);
    say('lessDistorted', verdict.lessDistorted, context);
    say('breathing', verdict.breathing, context);
    say('rigid', verdict.rigid, context);
    say('trigonalLean', verdict.trigonalLean, context);
    say('favouredAxes', verdict.favouredAxes, context);
    verdict.favouredAxes.families.forEach((family) => say('favouredAxes.families', family, {
        ...context,
        family: `<${family.family}>`,
        entry: data.local[`${family.axis}_axis_within_15deg`]?.[`<${family.family}>`],
        verb: AXIS_VERB[family.axis]
    }));
    say('splitMemberSites', verdict.splitMemberSites, context);
    say('window.bondsOutsideWindow', verdict.window.bondsOutsideWindow, context);
    say('window.contactsBetweenClusters', verdict.window.contactsBetweenClusters, context);
    // For the answer checks (review L5): the engine's answers.
    const facts = {
        cluster: {
            element,
            members: kind.members,
            window: [lo, hi],
            distorted: verdict.distorted.answer,
            static: verdict.staticDistortion.answer,
            breathing: verdict.breathing.answer,
            rigid: verdict.rigid.answer,
            split: verdict.splitMemberSites.answer,
            families: favoured.map((family) => `${AXIS_VERB[family.axis]} <${family.family}>`)
        },
        // For the single-configuration check (the user's rule, 2026-10-24):
        // the clusters' shape distribution in the run's one configuration.
        distribution: { kind: 'cluster_shape', subject: `the ${element}${kind.members} cluster-shape distribution`, configurations: 1 }
    };
    return { data, checks, facts };
};

// --- analyze_modulations ---------------------------------------------------------

// analyze_modulations' own result budget (characters): the rows of four
// elements with an ordered point's polarisation and their diffuse rows, and
// its checks. Since plan item 4.4's diffuse rows and note the largest results
// are three elements each: the Demo's 3,687 (three long-wave notes and the
// in-phase note), dev-zone-boundary's 3,344 (an ordered point); a fourth
// element adds about 500.
export const MODULATION_RESULT_BUDGET = 4500;

const qLabel = (q) => `(${q.map((value) => {
    const rounded = Math.round(value * 1000) / 1000;
    return rounded === 0.5 ? '½' : String(rounded + 0);
}).join(', ')})`;

// A mode row is [q label, ratio, percentile, sigma]; sigma reads '>8' where
// the engine reports σ = +∞ (as its rules' sigmaMax), and a ratio or σ that
// sits below its bar closer than six significant digits reads '<3' (F2).
const MODE_ROW = ['string', 'number|string|null', 'number|null', 'number|string|null'];
// An ordered point's displacement direction (plan item 4.3): the engine's
// polarisation as a [u v w] to two decimals and its share of the power;
// nearest_low_index only within LOW_INDEX_CONE_DEG of <100>, <110> or <111>.
const POLARISATION_ROW = {
    q: 'string',
    uvw: ['number', 'number', 'number'],
    share: 'number',
    'nearest_low_index?': { uvw: ['number', 'number', 'number'], off_deg: 'number|null' }
};
// The diffuse test (plan item 4.4) as a row: [lag, correlation, sigma] — the
// neighbour-cell lag whose correlation the power fits most strongly, that
// correlation, the whole neighbour-cell pattern's σ; absent when the grid
// has too few compared modes.
const DIFFUSE_ROW = ['string', 'number|string|null', 'number|string|null'];
const MODULATIONS_SHAPE = {
    note: 'string',
    supercell: ['number', 'number', 'number'],
    wavevectors: 'number',
    elements: [{
        element: 'string', sites: 'number', zone_boundary: [MODE_ROW], strongest: [MODE_ROW], 'polarisation?': [POLARISATION_ROW],
        'diffuse?': DIFFUSE_ROW
    }],
    'untestable_zone_boundary?': ['string']
};

// F2 (2026-10-19): a ratio or σ below its bar (ratio 3, 3σ: the order rule
// the note states) is never printed rounded onto it — "3σ" beside "≥ 3 is
// order" and a "no order" verdict. Where its usual rounding reads at or above
// the bar, it is shown with the fewest significant digits (3 up to
// BAR_DIGITS_MAX) that stay below it (2.96, not 3), and nearer than that as
// below the bar ('<3'). A value at or above its bar rounds to at or above it.
const BAR_DIGITS_MAX = 6;
const shownBelowBar = (value, bar, usual) => {
    if (!(Number.isFinite(value) && value < bar && Number.isFinite(usual) && usual >= bar)) return { shown: usual, near: false };
    for (let digits = 3; digits <= BAR_DIGITS_MAX; digits += 1) {
        const finer = round(value, digits);
        if (finer < bar) return { shown: finer, near: true };
    }
    return { shown: null, near: true };
};
// A row cell: the shown value, or '<bar' when no digits keep it below.
const barCell = (value, bar, usual) => {
    const { shown, near } = shownBelowBar(value, bar, usual);
    return near && shown == null ? `<${bar}` : shown;
};

// The ratio and σ in a row. A null σ in a sentence is σ = −∞ (an undefined
// one is said in words), written as the tool always wrote it.
const ratioCell = (ratio, rules) => barCell(ratio, rules.ratio, round(ratio, 2));
const sigmaCell = (sigma, rules, bar = rules.sigma) => (sigma == null ? null : (sigma >= rules.sigmaMax ? '>8' : barCell(sigma, bar, roundFixed(sigma, 1))));
// σ in a sentence: past 8 as '> 8σ'; null (σ = −∞: 1 − (1 − p)^F rounds to 1,
// so σ < −8.2) as '< −8σ' — until F1 no nearest miss was ever that far below
// the background, and the sentence printed 'nullσ'.
const sigmaText = (sigma) => {
    if (sigma == null) return '< −8σ';
    return sigma > 8 ? '> 8σ' : `${round(sigma, 2)}σ`;
};
const sameQ = (a, b) => a.length === b.length && a.every((value, axis) => value === b[axis]);
// A neighbour-cell lag [u v w] as written: "[1 -1 0]".
const lagLabel = (lag) => `[${lag.join(' ')}]`;
// A signed correlation to two significant digits, never printed onto its bar
// (F2's rule on |Γ| against rules.diffuseCorrelation): nearer than six
// digits, '<0.1' (or '>-0.1' for a negative one).
const correlationCell = (correlation, rules) => {
    if (correlation == null) return null;
    const cell = barCell(Math.abs(correlation), rules.diffuseCorrelation, round(Math.abs(correlation), 2));
    if (typeof cell === 'string') return correlation < 0 ? `>-${rules.diffuseCorrelation}` : cell;
    return correlation < 0 ? -cell : cell;
};
const signedText = (value) => (value < 0 ? `−${Math.abs(value)}` : `+${value}`);

// No element has an ordered zone-boundary point: one sentence for all of
// them, naming the nearest miss (the largest σ over every element tested). A
// ratio or σ just below its bar says so (F2): "2.96σ after the correction:
// below the 3σ bar".
const noCellDoubling = ({ best, rules }) => {
    if (!best) return { level: 'ok', text: 'No order that doubles the cell: no zone-boundary point stands out of the diffuse background.' };
    const { point } = best;
    const ratio = shownBelowBar(point.ratio, rules.ratio, round(point.ratio, 2));
    const sigma = point.because === 'undefined_statistic' || !(point.sigma <= 8)
        ? { near: false } : shownBelowBar(point.sigma, rules.sigma, round(point.sigma, 2));
    let significance = 'significance undefined';
    if (point.because !== 'undefined_statistic') significance = sigma.near ? `${sigma.shown ?? `< ${rules.sigma}`}σ` : sigmaText(point.sigma);
    const bars = [ratio.near ? `${rules.ratio}×` : null, sigma.near ? `${rules.sigma}σ` : null].filter(Boolean);
    const ratioText = ratio.near && ratio.shown == null ? `< ${rules.ratio}` : ratio.shown;
    return {
        level: 'ok',
        text: 'No order that doubles the cell: no zone-boundary point stands out of the diffuse background '
            + `(the strongest, ${best.element} at ${qLabel(point.q)}, is ${ratioText}× the median, `
            + `${significance} after the correction${bars.length ? `: below the ${bars.join(' and ')} bar${bars.length > 1 ? 's' : ''}` : ''}).`
    };
};

// An ordered point's polarisation as a row: its [u v w] to two decimals, its
// share, and the nearest low-index direction within LOW_INDEX_CONE_DEG.
const polarisationRow = (point, unitCell) => {
    const { direction, uvw, share } = point.polarisation;
    const family = lowIndexOf(direction, unitCell);
    return { q: qLabel(point.q), uvw: uvw.map((value) => roundFixed(value, 2)), share: roundFixed(share, 2), ...(family ? { nearest_low_index: family } : {}) };
};
// ", polarised along [0 0 1] (2.1° off; 99 % of the power at q)" (plan item 4.3).
const polarisedText = (row) => {
    if (!row) return '';
    const share = `${Math.round(row.share * 100)} % of the power at q`;
    return row.nearest_low_index
        ? `, polarised along [${row.nearest_low_index.uvw.join(' ')}] (${row.nearest_low_index.off_deg}° off; ${share})`
        : `, polarised along [${row.uvw.join(' ')}] (${share})`;
};

/**
 * How analyze_modulations says each element finding of the modulation
 * engine's verdict (workers/modulations.js, modulations.py): (finding,
 * because) → a check, or null. A cell doubling is said once per ordered
 * point, and its absence once for every element (noCellDoubling). Every
 * `because` the engine can give has an entry (copilotTools.test.js checks).
 */
// "Ga", "Ga and Se", "Ga, Se and Ta".
const listText = (names) => (names.length > 1 ? `${names.slice(0, -1).join(', ')} and ${names.at(-1)}` : names[0]);

/**
 * How the diffuse verdict (the modulation engine's diffuseOrder finding) is
 * read (the local-order interview, 2026-10-24): only a correlation in
 * antiphase that passes its bar is order between cells; one in phase is
 * correlated motion — neighbouring cells moving together — never order, and
 * correlation all the same. The tool's facts and arm A+
 * (validation/armA.js, SWEEP_RULES.diffuseOrder) read these.
 */
export const isDiffuseOrder = (finding) => finding?.because === 'correlated_antiphase';
export const isCorrelatedMotion = (finding) => finding?.because === 'correlated_in_phase';

export const MODULATION_SENTENCES = {
    cellDoubling: {
        ordered_point: ({ element, point, family, polarisation }) => ({
            level: 'ok',
            text: `Ordered ${element} displacements at q = ${qLabel(point.q)}: ${round(point.ratio, 3)}× the median power, `
                + `${sigmaText(point.sigma)} against independent cells (corrected over ${family} zone-boundary points) `
                + '— a superstructure doubling the cell along '
                + `${['a', 'b', 'c'].filter((_, axis) => point.q[axis]).join(', ')}${polarisedText(polarisation)}.`
        }),
        no_ordered_point: noCellDoubling,
        undefined_statistic: noCellDoubling,
        untestable: noCellDoubling,
        no_displacements: null
    },
    generalOrder: {
        ordered_mode: ({ element, mode, family, polarisation }) => ({
            level: 'ok',
            text: `Correlated ${element} displacements between cells at q = ${qLabel(mode.q)}: ${round(mode.ratio, 3)}× the median power, `
                + `${sigmaText(mode.sigma)} against independent cells (corrected over ${family} compared modes) — order at a `
                + `general wavevector (incommensurate, or short-range correlation between cells), not a cell doubling${polarisedText(polarisation)}.`
        }),
        no_ordered_mode: null,
        undefined_statistic: null,
        no_general_mode: null,
        no_displacements: null
    },
    // What a pattern is evidence of is inference, method as in the domain
    // rules (llm copilot/method.js; decided 2026-10-19): arm C's method 'off'
    // sends the long-wave and in-phase notes without their causal
    // parenthetical, and with their limit ("not a … superstructure").
    longWave: {
        long_wave_dominates: ({ element, mode }) => methodCheck('note',
            `${element}'s strongest mode is a long-wavelength wave, q = ${qLabel(mode.q)} at ${round(mode.ratio, 2)}× the median: `
                + 'box-scale displacement correlations',
            { method: ' (common in RMC)' },
            ', not a crystallographic superstructure.'),
        below_ratio: null,
        strongest_not_long_wave: null,
        undefined_statistic: null,
        no_displacements: null
    },
    // Plan item 4.4: correlation between neighbouring cells spread over many
    // q. In antiphase it is short-range order, said per element; in phase it
    // is neighbouring cells moving together, correlated motion and not order
    // (the local-order interview, 2026-10-24: only antiphase or general-q
    // order past its bar is order), said once for every element.
    diffuseOrder: {
        correlated_antiphase: ({ element, diffuse, family }) => ({
            level: 'ok',
            text: `Short-range order between cells: ${element} displacements in neighbouring cells ${lagLabel(diffuse.lag)} apart move in antiphase `
                + `(correlation ${signedText(round(diffuse.correlation, 2))}; the neighbour-cell pattern ${sigmaText(diffuse.sigma)} against independent `
                + `cells, corrected over ${family} elements), the power spread over many q toward the zone boundary with no single mode `
                + 'standing out — order that alternates between cells over a short range, not a long-range superstructure.'
        }),
        correlated_in_phase: ({ elements, strongest, diffuse }) => methodCheck('note',
            `${listText(elements)} displacements are correlated in phase between neighbouring cells (${elements.length > 1 ? `the strongest, ${strongest}, ` : ''}`
                + `${signedText(round(diffuse.correlation, 2))} between cells ${lagLabel(diffuse.lag)} apart, ${sigmaText(diffuse.sigma)} against `
                + 'independent cells): the power beyond the long waves leans toward the zone centre — neighbouring cells moving together '
                + '(correlated motion',
            { method: '; as RMC\'s long waves and acoustic-like thermal motion give' },
            '), not order.'),
        below_correlation: null,
        below_sigma: null,
        undefined_statistic: null,
        too_few_modes: null,
        no_displacements: null
    }
};

/**
 * The modulation result the model reads: the engine's points as compact rows,
 * each ordered point's polarisation, and its verdict said through
 * MODULATION_SENTENCES. The significance family is every element of the
 * configuration, whichever this call names (the engine's `family`, F3).
 * `unitCell` (rows a, b, c in Å) names a polarisation's nearest low-index
 * direction.
 */
export const summarizeModulations = (payload, { unitCell = null } = {}) => {
    const checks = [...parseWarningCheck(payload)];
    const tested = (payload.elements || []).filter((entry) => entry.sites);
    const { family, rules } = payload;
    const familyElements = family.elements ?? tested.map((entry) => entry.element);
    // The ordered points (a zone-boundary point, or a general mode) and their polarisation rows.
    const orderedPoints = (entry) => [
        ...(entry.zoneBoundary || []).filter((point) => point.ordered === true),
        ...(entry.strongest || []).filter((mode) => mode.kind === 'general' && mode.ordered === true)
    ].filter((point) => point.polarisation);
    const polarisations = new Map(tested.map((entry) => [entry.element, orderedPoints(entry).map((point) => polarisationRow(point, unitCell))]));
    const polarisationOf = (element, q) => polarisations.get(element).find((row) => row.q === qLabel(q)) ?? null;
    const anyPolarisation = [...polarisations.values()].some((rows) => rows.length);
    // The diffuse row of an element whose grid holds enough compared modes for the test.
    const diffuseRow = (entry) => (entry.diffuse?.df != null && entry.diffuse.lag
        ? [lagLabel(entry.diffuse.lag), correlationCell(entry.diffuse.correlation, rules), sigmaCell(entry.diffuse.sigma, rules, rules.diffuseSigma)]
        : null);
    const anyDiffuse = tested.some((entry) => diffuseRow(entry));
    const data = {
        note: 'zone_boundary rows are [q, ratio, percentile, sigma]: an element\'s displacement power at wavevector q (fractions '
            + 'of the cell\'s reciprocal axes) over its median across all q of the supercell grid (the diffuse background), the '
            + 'share of the compared modes (one per ±q pair, long waves with every |q_i| ≤ 0.25 left out) with less power, and '
            + `the calibrated one-sided significance against independent cells (each cell's own correlations kept), corrected over all ${family.zoneBoundaryPoints} zone-boundary `
            + `points of the configuration's ${familyElements.length} elements, whichever the call names (≥ 3 with ratio ≥ 3 is order). `
            + 'strongest rows are [q, ratio, percentile, sigma]: the sigma of a general q (neither zone boundary nor long wave; a '
            + 'complex sum, corrected over the compared modes of every element), null otherwise. '
            + 'A ½ component is order that doubles the cell along that axis; order at a general q is correlation between cells '
            + '(incommensurate or short-range) that does not.'
            + (anyDiffuse
                ? ' diffuse is [lag, correlation, sigma]: the strongest correlation between a site\'s displacements in cells one lag '
                    + 'apart that the power over the compared modes shows (+ in phase, − antiphase), and the significance of the 13 '
                    + 'neighbour-cell lags together against independent cells (≥ 3 with |correlation| ≥ 0.1: correlation spread over many q).'
                : '')
            + (anyPolarisation
                ? ' polarisation: each ordered point\'s displacement direction (the leading axis of the element\'s power tensor at '
                    + 'q) as a direct-lattice [u v w], an axis, and share, the fraction of the power at q along it (1: one '
                    + 'direction; ⅓: none).'
                : ''),
        supercell: payload.supercell,
        wavevectors: payload.wavevectors,
        elements: tested.map((entry) => ({
            element: entry.element,
            sites: entry.sites,
            // [q, ratio to the median, percentile, σ] — compact rows: 7 points × elements.
            zone_boundary: entry.zoneBoundary.map((point) => [
                qLabel(point.q), ratioCell(point.ratio, rules), roundFixed(point.percentile, 1), sigmaCell(point.sigma, rules)
            ]),
            strongest: entry.strongest.slice(0, 3).map((mode) => [
                qLabel(mode.q), ratioCell(mode.ratio, rules), roundFixed(mode.percentile, 1), sigmaCell(mode.sigma, rules)
            ]),
            ...(polarisations.get(entry.element).length ? { polarisation: polarisations.get(entry.element) } : {}),
            ...(diffuseRow(entry) ? { diffuse: diffuseRow(entry) } : {})
        }))
    };
    if (payload.untestable?.length) data.untestable_zone_boundary = payload.untestable.map((flags) => qLabel(flags.map((flag) => flag / 2)));
    if (!data.elements.length) {
        checks.push({ level: 'fail', text: 'No site has displacements to analyse (zero spread everywhere).' });
        return { data, checks };
    }
    const say = sayInto(checks, MODULATION_SENTENCES, 'analyze_modulations');
    const doubling = tested.filter((entry) => entry.verdict.cellDoubling.answer === true);
    if (doubling.length) {
        doubling.forEach((entry) => entry.verdict.cellDoubling.decidedBy.forEach((q) => say('cellDoubling', entry.verdict.cellDoubling, {
            element: entry.element,
            point: entry.zoneBoundary.find((point) => sameQ(point.q, q)),
            family: family.zoneBoundaryPoints,
            polarisation: polarisationOf(entry.element, q)
        })));
    } else {
        // Each element's nearest miss; the largest σ of them, the first element on a tie.
        const best = tested
            .map((entry) => ({ element: entry.element, point: entry.verdict.cellDoubling.nearest }))
            .filter((item) => item.point)
            .reduce((top, item) => (top == null || (item.point.sigma ?? -Infinity) > (top.point.sigma ?? -Infinity) ? item : top), null);
        say('cellDoubling', tested[0].verdict.cellDoubling, { best, rules });
    }
    // Correlation between cells at a general wavevector: each element's first ordered mode.
    tested.forEach((entry) => {
        const finding = entry.verdict.generalOrder;
        say('generalOrder', finding, {
            element: entry.element,
            mode: finding.answer === true ? entry.strongest.find((mode) => sameQ(mode.q, finding.decidedBy[0])) : null,
            family: family.comparedModes,
            polarisation: finding.answer === true ? polarisationOf(entry.element, finding.decidedBy[0]) : null
        });
    });
    tested.forEach((entry) => say('longWave', entry.verdict.longWave, { element: entry.element, mode: entry.strongest[0] }));
    // Diffuse order between cells (plan item 4.4): each element in antiphase
    // on its own; every element in phase in one note, naming the most
    // significant. An engine without the test (no diffuseOrder) says nothing.
    const diffuseBecause = (entry) => entry.verdict.diffuseOrder?.because;
    tested.filter((entry) => diffuseBecause(entry) === 'correlated_antiphase').forEach((entry) => say('diffuseOrder', entry.verdict.diffuseOrder, {
        element: entry.element, diffuse: entry.diffuse, family: family.diffuse
    }));
    const inPhase = tested.filter((entry) => diffuseBecause(entry) === 'correlated_in_phase');
    if (inPhase.length) {
        const top = inPhase.reduce((best, entry) => ((entry.diffuse.sigma ?? -Infinity) > (best.diffuse.sigma ?? -Infinity) ? entry : best));
        say('diffuseOrder', top.verdict.diffuseOrder, { elements: inPhase.map((entry) => entry.element), strongest: top.element, diffuse: top.diffuse });
    }
    tested.filter((entry) => entry.verdict.diffuseOrder && !['correlated_antiphase', 'correlated_in_phase'].includes(diffuseBecause(entry)))
        .forEach((entry) => say('diffuseOrder', entry.verdict.diffuseOrder, {}));
    if (payload.untestable?.length) {
        checks.push({ level: 'note', text: `An odd supercell axis cannot hold ${data.untestable_zone_boundary.join(', ')}: not tested.` });
    }
    // For the answer check: a superstructure claim needs an ordered point.
    // A claim of order at a box-scale wavevector (a long wave: every
    // |q_i| ≤ rules.longWaveQ, the engine's isLongWave) is flagged whatever
    // the mode's power (the user's rule, 2026-10-24): the rule's bound and the
    // long waves among each element's strongest modes, as the data shows them.
    const facts = {
        modulation: {
            elements: tested.map((entry) => entry.element),
            ordered: tested.flatMap((entry) => entry.verdict.cellDoubling.decidedBy.map((q) => ({ element: entry.element, q }))),
            general: tested
                .filter((entry) => entry.verdict.generalOrder.answer === true)
                .map((entry) => ({ element: entry.element, q: entry.verdict.generalOrder.decidedBy[0] })),
            long_wave_q: rules.longWaveQ,
            long_waves: tested.flatMap((entry) => entry.strongest.slice(0, 3)
                .filter((mode) => mode.kind === 'long_wave')
                .map((mode) => ({ element: entry.element, q: mode.q }))),
            // The diffuse verdicts that passed their bar (the local-order
            // interview, 2026-10-24): antiphase is order between cells, in
            // phase correlated motion — correlation, never order.
            diffuse: tested
                .filter((entry) => isDiffuseOrder(entry.verdict.diffuseOrder) || isCorrelatedMotion(entry.verdict.diffuseOrder))
                .map((entry) => ({
                    element: entry.element,
                    phase: isDiffuseOrder(entry.verdict.diffuseOrder) ? 'antiphase' : 'in_phase',
                    lag: [...entry.diffuse.lag],
                    correlation: round(entry.diffuse.correlation, 2),
                    sigma: round(entry.diffuse.sigma, 3)
                }))
        }
    };
    return { data, checks, facts };
};

// --- detect_symmetry -------------------------------------------------------------

// No average sites: tolerance_A only. A structure the finder did not analyse
// reads space_group 'not analysed' and n_ops '—'.
const SYMMETRY_SHAPE = {
    tolerance_A: 'number',
    'space_group?': 'string',
    'number?': 'number|null',
    'point_group?': 'string|null',
    'point_group_facts?': 'string',
    'n_ops?': 'number|string',
    'max_residual_A?': 'number|null',
    'orbits?': [{ element: 'string', wyckoff: 'string', site_symmetry: 'string', frac: oneOf(XYZ, 'null') }]
};

export const summarizeSymmetry = (described, tolerance) => {
    if (!described) {
        return {
            data: { tolerance_A: tolerance },
            checks: [{ level: 'fail', text: 'The structure has no average sites to analyse.' }]
        };
    }
    const facts = pointGroupFacts(described.pointGroup);
    const data = {
        tolerance_A: tolerance,
        space_group: described.spaceGroup,
        number: described.spaceGroupNumber ?? null,
        point_group: described.pointGroup,
        ...(facts ? { point_group_facts: describePointGroup(described.pointGroup) } : {}),
        n_ops: described.nSpace,
        max_residual_A: round(described.maxResidual, 2),
        orbits: (described.orbits || []).map((orbit) => ({
            element: orbit.element,
            wyckoff: orbitLabel(orbit),
            site_symmetry: orbit.site,
            frac: Array.isArray(orbit.rep) ? orbit.rep.map((value) => roundFixed(value, 4)) : null
        }))
    };
    const checks = [];
    if (described.skipped) {
        checks.push({ level: 'warn', text: `Not analysed: ${described.reason}` });
    } else if (described.spaceGroup === 'undetermined') {
        checks.push({ level: 'note', text: 'No closed group was found at this tolerance.' });
    } else if (String(described.spaceGroup).startsWith('≥')) {
        checks.push({ level: 'note', text: 'Only a lower bound (crystal class): no number, no Wyckoff letters.' });
    } else {
        checks.push({ level: 'ok', text: `${described.spaceGroup}${described.spaceGroupNumber ? ` (No. ${described.spaceGroupNumber})` : ''} holds at ${tolerance} Å.` });
    }
    // Where to read inversion, chirality and polarity from is method (llm
    // copilot/method.js): arm C's method 'off' sends the point group's facts alone.
    if (facts) {
        checks.push(methodCheck('note',
            `Point group ${describePointGroup(described.pointGroup)}.`,
            { method: ' Take inversion, chirality and polarity from this, not from memory.' }));
    }
    // For the answer check (llm/copilot/verification/facts.js): the group and what its
    // point group implies.
    // Its operations in the cell and its tolerance, for the group–subgroup
    // check (the local-symmetry interview, 2026-10-24).
    const result = { data, checks };
    if (facts && !String(described.spaceGroup).startsWith('≥') && described.spaceGroup !== 'undetermined') {
        result.facts = {
            groups: [{ symbol: described.spaceGroup, ...facts, ...(Number.isFinite(described.nSpace) ? { n_ops: described.nSpace } : {}), tolerance_A: tolerance }]
        };
    }
    return result;
};

// --- measure_symmetry_lowering ---------------------------------------------------

// Shifts reported per element: the largest site's, and how many sites move.
const SHIFT_ROWS = 12;
// A direction within this angle of a <100>, <110> or <111> direction is named
// by it: half the smallest angle between two such families (<110>–<111>,
// 35.3°), the validation's direction tolerance.
export const LOW_INDEX_CONE_DEG = 15;
// The nearest low-index direction of a Cartesian vector, within the cone, or null.
const lowIndexOf = (vector, unitCell) => {
    const found = unitCell ? nearestLowIndex(vector, unitCell, { maxIndex: 1, toleranceDeg: LOW_INDEX_CONE_DEG }) : null;
    return found?.exact ? { uvw: found.uvw, off_deg: roundFixed(found.offDeg, 1) } : null;
};

// {} when there is nothing to measure (no average sites, or not analysed);
// above_noise_floor only with a noise floor, nearest_low_index only for a
// shift within LOW_INDEX_CONE_DEG of <100>, <110> or <111>.
const LOWERING_SHAPE = oneOf({}, {
    note: 'string',
    parent: { tolerance_A: 'number', space_group: 'string', point_group: 'string|null' },
    lower: oneOf({ tolerance_A: 'number', space_group: 'string', point_group: 'string|null', 'point_group_facts?': 'string' }, 'null'),
    noise_floor_A: 'number|null',
    sites: [{
        ref: 'number',
        element: 'string',
        shift_A: 'number|null',
        direction_uvw: oneOf(XYZ, 'null'),
        'nearest_low_index?': { uvw: ['number', 'number', 'number'], off_deg: 'number|null' },
        se_A: 'number|null',
        'above_noise_floor?': 'boolean'
    }]
});

/**
 * The coherent shift of a symmetry lowering (symmetryShift.js), for the model:
 * the parent group and the lower group, each site's shift from the parent
 * positions with its crystal direction, against the sampling noise floor
 * (the review's L4, approved 2026-10-14). `lower` is describeSymmetry at
 * `tolerance`; `noiseFloor` the run context's (3√2 × rms/√copies).
 */
export const summarizeSymmetryLowering = (shifts, lower, { tolerance, noiseFloor = null, unitCell = null } = {}) => {
    if (!shifts) {
        return { data: {}, checks: [{ level: 'fail', text: 'The structure has no average sites to analyse.' }] };
    }
    if (shifts.skipped) return { data: {}, checks: [{ level: 'warn', text: `Not analysed: ${shifts.reason}` }] };
    const floor = Number.isFinite(noiseFloor) ? noiseFloor : null;
    const parentFacts = pointGroupFacts(shifts.parent.pointGroup);
    const lowerFacts = lower?.pointGroup ? pointGroupFacts(lower.pointGroup) : null;
    const data = {
        note: 'shift_A = the site mean\'s distance from the position the parent group puts it at, measured with the parent '
            + 'placed on the majority of the atoms (the origin a lowering leaves free is fixed by the geometric median of all '
            + 'shifts); a shift every copy of a site shares, which the displacement analyses cannot show; direction_uvw is a '
            + 'direct-lattice direction (±); se_A the sampling error of the site mean along it',
        parent: { tolerance_A: shifts.parent.tolerance, space_group: shifts.parent.spaceGroup, point_group: shifts.parent.pointGroup },
        lower: lower ? { tolerance_A: tolerance, space_group: lower.spaceGroup, point_group: lower.pointGroup } : null,
        noise_floor_A: floor == null ? null : round(floor, 2),
        sites: [...shifts.sites].sort((a, b) => b.length - a.length).slice(0, SHIFT_ROWS).map((site) => ({
            ref: site.referenceNumber,
            element: site.element,
            shift_A: roundFixed(site.length, 3),
            direction_uvw: site.length > 0 && unitCell ? cartesianToUvw(site.shift.map((v) => v / site.length), unitCell) : null,
            ...(site.length > 0 && lowIndexOf(site.shift, unitCell) ? { nearest_low_index: lowIndexOf(site.shift, unitCell) } : {}),
            se_A: site.standardError == null ? null : round(site.standardError, 2),
            ...(floor != null ? { above_noise_floor: site.length > floor } : {})
        }))
    };
    if (lowerFacts) data.lower.point_group_facts = describePointGroup(lower.pointGroup);
    const checks = [];
    if (!shifts.sites.length) {
        checks.push({ level: 'fail', text: `No operation of a parent group fits at ${shifts.parent.tolerance} Å.` });
        return { data, checks };
    }
    // The largest shift of each element, above the floor.
    const byElement = new Map();
    shifts.sites.forEach((site) => {
        const current = byElement.get(site.element);
        if (!current || site.length > current.length) byElement.set(site.element, site);
    });
    const moved = [...byElement.values()]
        .filter((site) => (floor == null ? false : site.length > floor))
        .sort((a, b) => b.length - a.length);
    const sameGroup = lower && lower.spaceGroup === shifts.parent.spaceGroup;
    const groups = lower && !sameGroup ? ` (the lower group ${lower.spaceGroup} at ${tolerance} Å` : '';
    if (moved.length) {
        moved.forEach((site) => {
            const uvw = unitCell ? cartesianToUvw(site.shift.map((v) => v / site.length), unitCell) : null;
            const family = lowIndexOf(site.shift, unitCell);
            const along = family ? ` along [${family.uvw.join(' ')}] (${family.off_deg}° off)` : (uvw ? ` along [${uvw.join(' ')}]` : '');
            const count = shifts.sites.filter((s) => s.element === site.element && floor != null && s.length > floor).length;
            checks.push({
                level: 'ok',
                text: `${site.element} shifted ${roundFixed(site.length, 3)} Å${along} from the ${shifts.parent.spaceGroup} `
                    + `positions (site ${site.referenceNumber}${count > 1 ? ` and ${count - 1} more` : ''}; noise floor ${round(floor, 2)} Å)${groups ? `${groups})` : ''}: `
                    + 'a coherent shift every copy shares, a static distortion of the average structure.'
            });
        });
        if (lowerFacts && parentFacts) {
            const changes = [
                parentFacts.centrosymmetric && !lowerFacts.centrosymmetric ? 'loses the inversion centre' : null,
                !parentFacts.polar && lowerFacts.polar ? 'becomes polar' : null
            ].filter(Boolean);
            if (changes.length && !sameGroup) {
                checks.push({ level: 'note', text: `${shifts.parent.spaceGroup} → ${lower.spaceGroup} ${changes.join(' and ')}.` });
            }
        }
    } else if (floor != null) {
        checks.push({
            level: 'ok',
            text: `No site is shifted from the ${shifts.parent.spaceGroup} positions above the ${round(floor, 2)} Å noise floor `
                + `(largest ${roundFixed(Math.max(...shifts.sites.map((s) => s.length)), 3)} Å): no coherent lowering of the average structure.`
        });
    } else {
        checks.push({ level: 'note', text: 'No noise floor for this structure: compare the shifts with their se_A.' });
    }
    if (sameGroup && moved.length) {
        checks.push({ level: 'note', text: `The lower tolerance still gives ${lower.spaceGroup}: the shifts lie within it.` });
    }
    const facts = {
        lowering: {
            parent: shifts.parent.spaceGroup,
            lower: lower?.spaceGroup ?? null,
            shifted: moved.map((site) => ({
                element: site.element,
                shift_A: roundFixed(site.length, 3),
                direction_uvw: lowIndexOf(site.shift, unitCell)?.uvw
                    ?? (unitCell ? cartesianToUvw(site.shift.map((v) => v / site.length), unitCell) : null)
            }))
        }
    };
    return { data, checks, facts };
};

// --- the tool list ---------------------------------------------------------------

const elementProperty = (elements, description) => (
    elements.length
        ? { type: 'string', enum: elements, description }
        : { type: 'string', description }
);

/**
 * The analysis tools for one run. `request(kind, params, { signal })` runs an
 * engine request (worker or Flask) and stops when `signal` aborts; null means
 * no engine is reachable (the static build with no run). `structure` is the
 * Dashboard's parsed structure; its average-structure `basis` (a browser-loaded
 * run's structureFromRmc6f, or a server folder's /api/structure; none for a
 * coords-only file) enables detect_symmetry and measure_symmetry_lowering, and
 * so does a server folder's `basisOmitted` (above the finder's site cap, where
 * both answer "not analysed", as for the browser's basis). `plotFiles` and
 * `runSettings` are the run's (CopilotChat's props of the same names): with
 * the structure they give the run's shell catalogue (llm/context/shells.js),
 * against which compute_bond_angles judges its windows by the rule the
 * Workbench asks on; without a partials file there is none.
 * Model arguments reach the engines only through the explicit
 * mapping in each run(), never by spreading.
 */
export const buildCopilotTools = ({ request, structure, elements = [], plotFiles = null, runSettings = null }) => {
    const tools = [];
    const unitCellFrom = (payload) => unitCellVectors(payload?.latticeVectors ?? structure?.latticeVectors, payload?.supercell ?? structure?.supercell);
    // The run's one shell reading, memoised per partials data, structure and
    // run settings: the catalogue the Workbench's questions read for the same run.
    const shells = () => shellCatalogue({ plotFiles, structure, runSettings });

    if (request) {
        tools.push({
            name: 'list_sites',
            description: 'List every reference site of the RMC configuration with its displacement statistics from a '
                + 'PCA of its displacement cloud: U_iso (Å²), principal RMS amplitudes (Å), anisotropy, Mardia and '
                + 'per-axis kurtosis, PC1 as [u v w], and a split site\'s static displacement (Å) with its sampling error.',
            // The method advice (llm copilot/method.js): sent with method 'full', after the description.
            method: 'Use it to find the most disordered, anisotropic or anharmonic sites (their `ref`).',
            parameters: {
                type: 'object',
                properties: {
                    element: elementProperty(elements, 'Only the sites holding this element.')
                }
            },
            dataShape: SITES_SHAPE,
            label: ({ element }) => `Site displacement table${element ? ` · ${element}` : ''}`,
            page: 'ellipsoids',
            cards: () => ['pca-stats'],
            run: async ({ element }, { signal } = {}) => {
                const payload = await request('sites', {}, { signal });
                return summarizeSites(payload, { element, unitCell: unitCellFrom(payload) });
            }
        });

        tools.push({
            name: 'analyze_displacement_directions',
            description: 'Histogram the displacement directions of one site (site_ref, from list_sites) or of every atom '
                + 'of one element (element; "all" pools every atom) on the sphere, and test them: the peak direction '
                + '(Cartesian and crystal [u v w], ± its cell\'s radius) with its calibrated significance, labelled a finding, a hint '
                + 'or nothing (a finding is checked at a neighbouring resolution), the whole-map test against '
                + 'isotropy, the +u/−u asymmetry, and the orientation-tensor anisotropy.',
            parameters: {
                type: 'object',
                properties: {
                    site_ref: { type: 'integer', minimum: 1, description: 'Reference site number.' },
                    element: elementProperty(elements.length ? [...elements, 'all'] : [], 'Pool every atom of this element instead of one site.'),
                    frequency: { type: 'integer', minimum: 1, maximum: 24, description: 'Sphere resolution (10ν²+2 cells). Default: chosen from the point count.' }
                }
            },
            dataShape: ORIENTATION_SHAPE,
            label: ({ site_ref: ref, element }) => `Displacement directions ${ref != null ? `site ${ref}` : (element || '')}`.trim(),
            page: 'orientation',
            cards: () => ['direction-sphere'],
            // The page shows one site at a time (no element pooling); no
            // frequency means the engine chose it, the page's 'auto'.
            pageParams: ({ site_ref: ref, frequency }) => (ref != null
                ? { referenceNumber: ref, frequency: frequency ?? 'auto' }
                : null),
            run: async ({ site_ref: referenceNumber, element, frequency }, { signal } = {}) => {
                if (referenceNumber == null && !element) throw new Error('give site_ref (from list_sites) or element');
                const params = { geometry: false };
                if (referenceNumber != null) params.referenceNumber = referenceNumber;
                else params.element = element;
                if (frequency != null) params.frequency = frequency;
                const payload = await request('orientation', params, { signal });
                // A finding is judged at the engine's own resolution and holds
                // at the neighbouring one (the user's tiers, 2026-10-24): a
                // second engine call, only for a candidate (≥ 3σ, ≥ 1.2×).
                const auto = frequency == null || frequency === payload.recommendedFrequency;
                const neighbour = auto && directionCandidate(payload)
                    ? await request('orientation', { ...params, frequency: neighbourFrequency(payload.frequency) }, { signal })
                    : null;
                const basis = structure?.basis || [];
                const pooledSites = referenceNumber != null
                    ? 0
                    : (element === 'all' ? basis.length : basis.filter((site) => site.el === element).length);
                return summarizeOrientation(payload, unitCellFrom(payload), { pooledSites, neighbour, auto });
            }
        });

        tools.push({
            name: 'compute_bond_angles',
            description: 'Bond-angle distribution A–B–C with B the central atom (RMCProfile triplets): every pair of '
                + 'B neighbours, A within [r12_min, r12_max] Å of B and C within [r23_min, r23_max] Å (default: the '
                + 'A–B window). Returns the angle modes, mean ± sd, a 10° histogram, bond lengths and the B '
                + 'coordination, and for about six or four neighbours the split against an ideal octahedron or '
                + 'tetrahedron.',
            // The refusal itself stays in window_confirmed's own description.
            method: 'Use the whole first coordination shell (find_gr_peaks\' first_shell.window_A): a window '
                + 'that cuts the shell is put to the user first, and refused unless the user chose it (window_confirmed).',
            parameters: {
                type: 'object',
                properties: {
                    end1: elementProperty(elements, 'Element A.'),
                    apex: elementProperty(elements, 'Element B, the central atom.'),
                    end2: elementProperty(elements, 'Element C.'),
                    r12_min: { type: 'number', minimum: 0, maximum: 15, description: 'A–B window minimum (Å).' },
                    r12_max: { type: 'number', minimum: 0, maximum: 15, description: 'A–B window maximum (Å).' },
                    r23_min: { type: 'number', minimum: 0, maximum: 15, description: 'B–C window minimum (Å); give both or neither.' },
                    r23_max: { type: 'number', minimum: 0, maximum: 15, description: 'B–C window maximum (Å).' },
                    bin_width: { type: 'number', minimum: 0.05, maximum: 30, description: 'Histogram bin (degrees). Default 1.' },
                    window_confirmed: {
                        type: 'boolean',
                        description: 'true only when the user chose or typed this window (in the question, or through ask_user). '
                            + 'A window that cuts the coordination shell is refused unless it is set.'
                    }
                },
                required: ['end1', 'apex', 'end2', 'r12_min', 'r12_max']
            },
            // The bond windows, for the Workbench's check that a window does
            // not cut a coordination shell before this runs (llm clarify.js),
            // and for this tool's own refusal: the same windows, the same rule.
            shellWindows: bondShellWindows,
            dataShape: BOND_ANGLES_SHAPE,
            label: ({ end1, apex, end2, r12_min: lo, r12_max: hi }) => (
                `Bond angles ${end1}–${apex}–${end2}${Number.isFinite(lo) && Number.isFinite(hi) ? ` · ${lo}–${hi} Å` : ''}`
            ),
            page: 'geometry',
            cards: () => ['bond-angles'],
            pageParams: (args) => ({
                end1: args.end1,
                apex: args.apex,
                end2: args.end2,
                r12Min: args.r12_min,
                r12Max: args.r12_max,
                r23Min: args.r23_min ?? null,
                r23Max: args.r23_max ?? null,
                binWidth: args.bin_width ?? 1
            }),
            run: async (args, { signal } = {}) => {
                if ((args.r23_min == null) !== (args.r23_max == null)) {
                    throw new Error('give r23_min and r23_max together, or neither');
                }
                if (args.r12_min >= args.r12_max) throw new Error('r12_min must be below r12_max');
                if (args.r23_min != null && args.r23_min >= args.r23_max) throw new Error('r23_min must be below r23_max');
                const params = {
                    end1: args.end1,
                    apex: args.apex,
                    end2: args.end2,
                    r12Min: args.r12_min,
                    r12Max: args.r12_max,
                    binWidth: args.bin_width ?? 1
                };
                if (args.r23_min != null) {
                    params.r23Min = args.r23_min;
                    params.r23Max = args.r23_max;
                }
                const payload = await request('triplets', params, { signal });
                // Each window the Workbench reads, against the run's shell of its pair.
                const catalogue = shells();
                const readings = bondShellWindows(args).map(({ centre, partner, window }) => catalogue?.window(centre, partner, window) ?? null);
                return summarizeBondAngles(payload, { windowConfirmed: args.window_confirmed === true, shells: readings });
            }
        });

        tools.push({
            name: 'analyze_clusters',
            description: 'Shape of the finite clusters one element forms — the M4 tetrahedra of a lacunar spinel, '
                + 'dimers, trimers — with bonds in [r_min, r_max] Å in the average structure. Every cluster in the '
                + 'configuration is compared with a shuffle null in which each atom keeps its own site\'s displacements but takes them from '
                + 'a random other cluster: is there a distortion beyond independent atoms, which edge pattern '
                + '(tetrahedra: E/T2 modes, trigonal vs one edge), along which lattice axis (<100>/<110>/<111>), and is '
                + 'the average cluster itself distorted.',
            // Where the window comes from (inside the first sentence), and which
            // questions need the tool: set below, once the list shows whether
            // the run offers polyhedral_distortion (clusterMethod).
            method: clusterMethod({ polyhedra: false }),
            parameters: {
                type: 'object',
                properties: {
                    element: elementProperty(elements, 'Element whose atoms form the clusters.'),
                    r_min: { type: 'number', minimum: 0, maximum: 8, description: 'Bond window minimum (Å).' },
                    r_max: { type: 'number', minimum: 0.5, maximum: 8, description: 'Bond window maximum (Å).' }
                },
                required: ['element', 'r_min', 'r_max']
            },
            shellWindows: (args) => [{
                centre: args.element,
                partner: args.element,
                window: [args.r_min, args.r_max],
                apply: ([lo, hi]) => ({ ...args, r_min: lo, r_max: hi })
            }],
            dataShape: CLUSTERS_SHAPE,
            label: ({ element, r_min: lo, r_max: hi }) => (
                `Clusters ${element || ''}${Number.isFinite(lo) && Number.isFinite(hi) ? ` · ${lo}–${hi} Å` : ''}`.trim()
            ),
            // The shape and size rows with their checks run to ~3.4k characters
            // (the demo's Ta4): within the default 3k they reached the model as
            // a cut-off text, whose numbers no reference could name (THY-18).
            resultBudget: CLUSTER_RESULT_BUDGET,
            run: async ({ element, r_min: rMin, r_max: rMax }, { signal } = {}) => {
                if (rMin >= rMax) throw new Error('r_min must be below r_max');
                // The engine's verdict holds the split member sites too: it
                // reads them from the clouds it holds.
                return summarizeClusters(await request('clusters', { element, rMin, rMax }, { signal }));
            }
        });

        tools.push({
            name: 'analyze_modulations',
            description: 'Look for order that enlarges the unit cell — a superstructure the folded average structure '
                + 'cannot show (cell doubling as in low-temperature phases, octahedral tilts, '
                + 'charge-density waves). Fourier-transforms every site\'s displacements over the supercell\'s cell grid and '
                + 'compares the zone-boundary points (½,0,0)…(½,½,½) and the strongest modes with the configuration\'s own '
                + 'diffuse background (ratio to the median power, percentile). Optionally one element.',
            // Which other tool cannot show it.
            method: [{ after: 'the folded average structure', text: ' and detect_symmetry' }],
            parameters: {
                type: 'object',
                properties: {
                    element: elementProperty(elements, 'Only this element (default: each element).')
                }
            },
            dataShape: MODULATIONS_SHAPE,
            label: ({ element }) => `Superstructure test${element ? ` · ${element}` : ''}`,
            // Three elements with an ordered point, its polarisation and the
            // diffuse rows run to ~3.3k characters with the checks
            // (dev-zone-boundary), a fourth element to ~3.9k: within the
            // default budget the rows would be cut.
            resultBudget: MODULATION_RESULT_BUDGET,
            run: async ({ element }, { signal } = {}) => summarizeModulations(
                await request('modulations', element ? { element } : {}, { signal }),
                { unitCell: unitCellFrom(null) }
            )
        });
    }

    // The browser offers both symmetry tools for any basis, and above the
    // finder's site cap they answer "not analysed". A server folder above the
    // cap sends only the basis's size (`basisOmitted`): the same two tools,
    // with the same answers.
    if (structure?.basis?.length || structure?.basisOmitted) {
        tools.push({
            name: 'detect_symmetry',
            description: 'Run the space-group finder on the average structure (the folded reference sites) at a '
                + 'position tolerance in Å, as the Detected SG card does: group, number (only when it is certain), '
                + 'operations, and the Wyckoff orbits. Does not change the card.',
            method: [{ after: 'and the Wyckoff orbits.', text: ' Compare tolerances to see which symmetry the local distortion breaks.' }],
            parameters: {
                type: 'object',
                properties: {
                    tolerance_A: { type: 'number', minimum: 0.001, maximum: 1, description: 'Position tolerance (Å); the card default is 0.2.' }
                },
                required: ['tolerance_A']
            },
            dataShape: SYMMETRY_SHAPE,
            label: ({ tolerance_A: tol }) => `Symmetry at ${tol} Å`,
            page: 'dashboard',
            cards: () => ['detected-sg'],
            run: async ({ tolerance_A: tol }) => summarizeSymmetry(describeSymmetry(structure, tol), tol)
        });

        // Without rows (`basisOmitted`) the noise floor is unknown: null, and
        // the default below assumes 0.02 Å, as for a basis whose sites have one
        // copy each (the usual case above the cap: a glass or a 1×1×1 box).
        const floor = symmetryNoiseFloor(structure.basis);
        // The lower group's default tolerance: NOISE_FLOOR_MARGIN × the noise
        // floor, just above it (the one rule; arm A+ reads it too).
        const lowerDefault = lowerGroupTolerance(structure.basis);
        tools.push({
            name: 'measure_symmetry_lowering',
            description: 'How far each average site sits from the position a higher-symmetry (parent) group puts it at: the '
                + 'size and direction of a coherent shift every copy of a site shares (a polar off-centring of the average '
                + 'structure). The parent is the group at parent_tolerance_A (default 0.2 Å, the Detected SG card\'s); '
                + 'the lower group is the one '
                + `at tolerance_A (default ${lowerDefault} Å, just above the sampling noise floor). Each shift is in Å with its `
                + 'crystal [u v w], measured with the parent on the majority of the atoms, against the noise floor.',
            // Which analyses cannot show it, and why.
            method: [{
                after: '(a polar off-centring of the average structure)',
                text: ', which the displacement analyses cannot show — each site\'s cloud is centred on its own mean'
            }],
            parameters: {
                type: 'object',
                properties: {
                    parent_tolerance_A: { type: 'number', minimum: 0.001, maximum: 1, description: 'Tolerance of the parent group (Å). Default 0.2.' },
                    tolerance_A: { type: 'number', minimum: 0.001, maximum: 1, description: `Tolerance of the lower group (Å). Default ${lowerDefault}.` }
                }
            },
            dataShape: LOWERING_SHAPE,
            label: ({ parent_tolerance_A: parentTol = 0.2 }) => `Symmetry lowering from ${parentTol} Å`,
            page: 'dashboard',
            cards: () => ['detected-sg'],
            run: async ({ parent_tolerance_A: parentTol = 0.2, tolerance_A: tol = lowerDefault }) => summarizeSymmetryLowering(
                symmetryShifts(structure, parentTol),
                describeSymmetry(structure, tol),
                { tolerance: tol, noiseFloor: floor, unitCell: unitCellVectors(structure.latticeVectors, structure.supercell) }
            )
        });
    }
    // analyze_clusters' method names polyhedral_distortion for an extended
    // network only where the run offers it (that tool is its own PR: with or
    // without it, the advice names no tool the model lacks).
    const clusters = tools.find((tool) => tool.name === 'analyze_clusters');
    if (clusters) clusters.method = clusterMethod({ polyhedra: tools.some((tool) => tool.name === 'polyhedral_distortion') });
    return tools;
};
