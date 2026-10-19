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
// The summarizers are pure and exported for tests.

import { describeSymmetry, orbitLabel } from './symmetryModel';
import { siteLabel } from './siteLabel';
import { crystalOrientationRows, invert3, matVec3, transpose3, unitCellVectors } from './pcaCrystalFrame';
import { describePointGroup, pointGroupFacts } from './llm';
import { symmetryNoiseFloor } from './llm/context/runContext';
import { symmetryShifts } from './symmetryShift';
import { nearestLowIndex } from './viewer3d';
import { normalQuantile, regularizedGamma } from './workers/orientation';
import { oneOf } from './llm/copilot/dataShape';

// Each tool's `dataShape` (llm/copilot/dataShape.js) is the shape of the data
// its summarizer builds — the paths the model quotes by reference — declared
// beside that summarizer and carried on the tool's spec.
const XYZ = ['number|null', 'number|null', 'number|null'];

// A tolerance this factor above the sampling noise floor holds clearly above
// it (the Workbench's tolerance question uses the same margin, clarify.js).
export const NOISE_FLOOR_MARGIN = 1.25;

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
// Central atoms with no A–B bond above this fraction: the window misses them.
export const UNBONDED_APEX_FRACTION = 0.5;
// Fewer central atoms than this share with the most common A–B bond count:
// the coordination is spread, as when the window edge cuts through a shell
// (demo Se–Ta–Se at 2.3–2.7 Å: 47 %; whole shells there: 83–96 %).
export const MODAL_COORDINATION_SHARE = 0.6;

const parseWarningCheck = (payload) => (
    payload?.parseWarning ? [{ level: 'warn', text: `Structure file: ${payload.parseWarning}` }] : []
);

// --- split sites (list_sites and the cluster summary) ------------------------------

// A site is split (two positions) when its PC1 is resolved and its PC1 excess
// kurtosis is at most −0.5 — the value at which two equal wells ±δ of thermal
// width σ separate (δ = σ: κ = −2δ⁴/(δ² + σ²)² = −0.5) — and at least this
// many sampling errors (√(24/n)) below zero.
export const SPLIT_SITE_KURTOSIS = -0.5;
export const SPLIT_SITE_SE = 3;
// The threshold as the texts the model reads write it (a true minus sign).
const SPLIT_RULE_KURTOSIS = String(SPLIT_SITE_KURTOSIS).replace('-', '−');

/** Whether one site (an engine 'sites' row) is split: the app's one rule, for list_sites and analyze_clusters. */
export const isSplitSite = (site) => {
    const kappa = site?.excessKurtosis?.[0];
    return Boolean(site?.axisResolved?.[0]) && Number.isFinite(kappa) && site.count > 0
        && kappa <= SPLIT_SITE_KURTOSIS && kappa <= -SPLIT_SITE_SE * Math.sqrt(24 / site.count);
};

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

/** The sites (list_sites' engine rows) among `refs` whose clouds are split. */
export const splitMemberSites = (siteRows, refs) => (siteRows || []).filter((site) => (
    refs.includes(site.referenceNumber) && isSplitSite(site)
));

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
    return { data, checks };
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

// direction_uvw, angles_to_abc_deg and principal_axis_uvw need the unit cell.
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
        'angles_to_abc_deg?': XYZ
    },
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

export const summarizeOrientation = (payload, unitCell, { pooledSites = 0 } = {}) => {
    const direction = payload.peakDirection;
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
            + 'the map is not antipodally folded, so +u and −u differ; [u v w] is a direct-lattice direction, not a plane normal',
        site: payload.referenceNumber ?? null,
        element: payload.element ?? 'all',
        frequency: payload.frequency,
        cells: payload.cellCount,
        used_points: payload.usedPoints,
        rejected_points: payload.rejectedPoints,
        enhancement_range: [round(payload.vmin), round(payload.vmax)],
        peak,
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
    if (payload.peakSignificance >= SIGNIFICANCE_SIGMA) {
        checks.push({
            level: 'ok',
            text: `Preferred direction: the peak cell is ${round(payload.peakEnhancement)}× isotropic at ${round(payload.peakSignificance)}σ (corrected over all cells).`
        });
    } else {
        checks.push({ level: 'note', text: `No significant peak direction (${round(payload.peakSignificance)}σ < ${SIGNIFICANCE_SIGMA}σ).` });
    }
    if (payload.mapSignificance == null) {
        checks.push({ level: 'note', text: 'The whole-map test is undefined here (too few points per cell).' });
    } else if (payload.mapSignificance >= SIGNIFICANCE_SIGMA) {
        // Any anisotropic Gaussian ellipsoid passes the map test (a
        // thermal-only O flattened across its bond: map 7.2σ, peak 2.3σ):
        // it shows anisotropy, and only the peak test a preferred direction.
        checks.push({
            level: 'ok',
            text: `The displacements are anisotropic: the map departs from isotropy at ${round(payload.mapSignificance)}σ. `
                + (payload.peakSignificance >= SIGNIFICANCE_SIGMA
                    ? 'The peak test above names the preferred direction.'
                    : 'Thermal ellipsoids pass this test too: it is not a preferred direction, which needs the peak test.')
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
    if (payload.referenceNumber == null && pooledSites > 1) {
        checks.push({
            level: 'note',
            text: `Pooled over ${pooledSites} sites. Where their local axes differ (the four vertices of a tetrahedral cluster), `
                + 'a direction each site prefers in its own frame cancels here, so an isotropic pooled map is not evidence '
                + 'against a per-site or cluster distortion. Test single sites (site_ref), or analyze_clusters for a cluster\'s shape.'
        });
    }
    // For the answer check: no "preferred direction" unless the peak test of
    // one of the question's analyses found it at ≥ 3σ. The map test alone is
    // anisotropy (review L3, approved 2026-10-14): before, it counted too, and
    // let a preferred direction through on plain thermal anisotropy.
    const facts = {
        direction: {
            subject: payload.referenceNumber != null ? `site ${payload.referenceNumber}` : (payload.element || 'all atoms'),
            significant: payload.peakSignificance >= SIGNIFICANCE_SIGMA,
            anisotropic: payload.mapSignificance >= SIGNIFICANCE_SIGMA,
            peak_sigma: round(payload.peakSignificance),
            map_sigma: payload.mapSignificance == null ? null : round(payload.mapSignificance)
        }
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
    if (Math.abs(meanCoordination - 6) <= 0.5) {
        return {
            ideal: 'octahedron: 80 % cis at 90°, 20 % trans at 180°',
            cis_below_135_deg: band((center) => center < POLYHEDRON_SPLIT_DEG),
            trans_above_135_deg: band((center) => center >= POLYHEDRON_SPLIT_DEG)
        };
    }
    if (Math.abs(meanCoordination - 4) <= 0.5) {
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
// A refused window keeps the counts that show the cut and drops the angles
// (mean_deg … histogram_10deg_sin_corrected), with `refused` saying why; the
// B–C window and bond read 'same as A–B' for one shared window; `polyhedron`
// only for a whole shell of about six (octahedron) or four (tetrahedron).
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
    coordination: { note: 'string', counts: ['number'], mean: 'number|null' },
    'refused?': 'string',
    'polyhedron?': oneOf(
        { ideal: 'string', cis_below_135_deg: ANGLE_BAND, trans_above_135_deg: ANGLE_BAND },
        { ideal: 'string', within_15_deg_of_109_5: ANGLE_BAND }
    )
};

/**
 * The bond-angle result the model reads, with deterministic checks: a window
 * that misses most central atoms, or in which fewer than 60 % of them share
 * the most common bond count, cuts the coordination shell (the 2026-10-10 spot
 * check read such a cut as "incomplete octahedra"). Such a window is refused —
 * its angles are not the polyhedron's — unless the user chose it
 * (`windowConfirmed`, the coordination-shells skill's firm rule; its 60 %
 * threshold is an interview item): the counts that show the cut are returned
 * instead. For a whole shell of about six or four neighbours the angles are
 * also read against an ideal octahedron or tetrahedron, and the sin
 * θ-corrected histogram shows the trans angles raw counts hide.
 */
export const summarizeBondAngles = (payload, { windowConfirmed = false } = {}) => {
    const total = payload.angleCount || 0;
    const coordination = Array.isArray(payload.coordination) ? payload.coordination : [];
    const apexTotal = coordination.reduce((sum, count) => sum + count, 0);
    const meanCoordination = apexTotal
        ? coordination.reduce((sum, count, bonds) => sum + count * bonds, 0) / apexTotal
        : null;
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
            mean: round(meanCoordination)
        }
    };
    const checks = [...parseWarningCheck(payload)];
    if (!total) {
        checks.push({
            level: 'fail',
            text: 'No angles: no central atom has two neighbours in the window(s). Widen the window '
                + '(find_gr_peaks gives a first-peak window) or check the elements.'
        });
        return { data, checks };
    }
    let cutText = null;
    const unbonded = apexTotal ? (coordination[0] || 0) / apexTotal : 0;
    if (unbonded > UNBONDED_APEX_FRACTION) {
        cutText = `${Math.round(unbonded * 100)} % of the ${b} atoms have no ${a}–${b} bond in the window: it misses most of the shell.`;
    } else if (apexTotal) {
        const modalCount = Math.max(...coordination);
        const modal = coordination.indexOf(modalCount);
        const share = modalCount / apexTotal;
        const present = coordination.map((count, bonds) => (count > 0 ? bonds : null)).filter((bonds) => bonds != null);
        if (share < MODAL_COORDINATION_SHARE) {
            cutText = `Only ${Math.round(share * 100)} % of the ${b} atoms have the most common count (${modal} ${a}–${b} bonds); `
                + `the counts spread from ${present[0]} to ${present.at(-1)}. The window cuts through a coordination shell.`;
        }
    }
    // For the coordination answer check (review L5): the A–B count in this
    // window, and whether it is the whole shell or a window the user chose.
    const facts = {
        coordination: {
            centre: b,
            partner: a,
            window: payload.bond12 ?? null,
            cn: round(meanCoordination),
            whole_shell: !cutText,
            confirmed: windowConfirmed
        }
    };
    if (cutText && !windowConfirmed) {
        // Refused: keep what shows the cut, never the angles.
        ['mean_deg', 'sd_deg', 'modes', 'histogram_10deg', 'histogram_10deg_sin_corrected'].forEach((key) => { delete data[key]; });
        data.refused = 'the window cuts the coordination shell, so its angles are not the polyhedron\'s and are not returned';
        checks.push({
            level: 'fail',
            text: `Refused: ${cutText} Use the whole shell (find_gr_peaks' first_shell.window_A). If the user wants this window, `
                + 'ask them (ask_user) and call again with window_confirmed: true.'
        });
        return { data, checks, facts };
    }
    if (cutText) checks.push({ level: 'warn', text: `${cutText} Kept: the user chose this window; say that its angles are of a partial shell.` });
    const cut = Boolean(cutText);
    if (!cut) {
        checks.push({ level: 'ok', text: `${total} angles around ${apexTotal} ${b} atoms (mean coordination ${round(meanCoordination)}).` });
    }
    // Only a whole shell is read against an ideal polyhedron: a cut TaSe6
    // shell holds ~3.9 Se per Ta, which is not a tetrahedron. And only one
    // shell: an A–B–A triplet with one window. A mixed (Se–Ta–Te) or split
    // (two windows) triplet's angles are a subset of the polyhedron's.
    const oneShell = a === c && payload.sharedEnds === true;
    const reference = cut || !oneShell ? null : polyhedronReference(payload.binCenters || [], payload.counts || [], total, meanCoordination);
    if (reference) data.polyhedron = reference;
    if (reference?.trans_above_135_deg) {
        // A band with no angles has no mean: say "none", never "null°".
        const meanOf = (band) => (band.mean_deg == null ? 'none' : `${band.mean_deg}°`);
        checks.push({
            level: 'note',
            text: `Raw counts deplete angles near 180° (as sin θ): ${Math.round(reference.trans_above_135_deg.share * 100)} % of these angles `
                + `lie above 135° (mean ${meanOf(reference.trans_above_135_deg)}) and ${Math.round(reference.cis_below_135_deg.share * 100)} % `
                + `below (mean ${meanOf(reference.cis_below_135_deg)}), against an ideal octahedron's 20 % trans at 180° and 80 % cis `
                + 'at 90°. Read the shape from histogram_10deg_sin_corrected (1 = random directions), not from the raw modes.'
        });
    }
    return { data, checks, facts };
};

// --- analyze_clusters ------------------------------------------------------------

// A shape statistic departs from independent atoms only when both its z and
// its size clear these: with thousands of clusters a 1 % effect reaches
// z > 3 and still means nothing physically.
export const CLUSTER_Z = 3;
// Relative change (observed / null − 1) of a width or an rms.
export const CLUSTER_EFFECT = 0.05;
// Change of the tetrahedral T2 trigonality (its range is 1/3 to 1).
export const TRIGONALITY_EFFECT = 0.02;
// Observed / null fraction of clusters with an axis near one lattice family.
export const AXIS_ENHANCEMENT = 1.1;
// The average cluster counts as statically distorted when its edges spread
// over more than this many sampling errors (six independent normal edges span
// ~2.5 on average and exceed 5 in under 1 % of cases).
export const STATIC_SPREAD_SE = 5;
// Template bonds outside the window, or contacts linking clusters, per atom.
export const CLUSTER_OUTLIER_FRACTION = 0.05;
// analyze_clusters' own result budget (characters of the tool message).
export const CLUSTER_RESULT_BUDGET = 4500;

const significant = (z, change, effect) => (
    Number.isFinite(z) && Math.abs(z) >= CLUSTER_Z && Number.isFinite(change) && Math.abs(change) >= effect
);

// A member site is split by the rule list_sites reports (isSplitSite, above).

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

/**
 * The cluster result the model reads. The shape statistics (edge spread,
 * shape rms, tetrahedral modes, gyration-axis families) are read against the
 * engine's shape null — independent atoms at each cluster's own breathing —
 * and the size statistics (mean edge, breathing) against independent atoms:
 * against the latter alone, a collective size fluctuation read as "LESS
 * distorted" and as a shape preference (the review's L2, approved
 * 2026-10-14). `siteRows` (the engine's 'sites' rows, when the tool could
 * read them) tell split member sites, for which the null already holds the
 * distortion's amplitude.
 */
export const summarizeClusters = (payload, { siteRows = null } = {}) => {
    const [lo, hi] = payload.window || [];
    const element = payload.element;
    const data = {
        element,
        window_A: [lo, hi],
        sites: payload.sites?.length ?? 0
    };
    if (payload.isolatedSites) data.sites_without_bond = payload.isolatedSites;
    const checks = [...parseWarningCheck(payload)];
    if (payload.extended) {
        checks.push({
            level: 'fail',
            text: `The ${lo}–${hi} Å window links the ${element} atoms into an extended network (chains or a 3-D net), `
                + `not finite clusters: narrow it to the intra-cluster bond (find_gr_peaks: first_shell.window_A of ${element}-${element}).`
        });
        return { data, checks };
    }
    if (!payload.templates?.length || !payload.kind) {
        checks.push({
            level: 'fail',
            text: `No two ${element} sites of the average structure lie ${lo}–${hi} Å apart: no clusters. `
                + `Widen the window (find_gr_peaks gives the ${element}-${element} first shell).`
        });
        return { data, checks };
    }
    const { observed: obs } = payload;
    // Shape statistics against the shape null (when the engine gave one).
    const shape = payload.shapeNull ? { mean: payload.shapeNull.mean, z: payload.shapeNull.z } : { mean: payload.nullMean, z: payload.z };
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
    const families = (name) => (obs[name] ? Object.fromEntries(Object.keys(obs[name]).map((family) => [
        `<${family}>`, row(obs[name][family], shape.mean[name][family], shape.z[name][family], 3)
    ])) : undefined);
    if (obs.longAxis) data.local.long_axis_within_15deg = families('longAxis');
    if (obs.shortAxis) data.local.short_axis_within_15deg = families('shortAxis');
    data.outliers = {
        bonds_outside_window: round(payload.outliers.brokenEdgeFraction, 2),
        contacts_between_clusters_per_atom: round(payload.outliers.extraContactsPerAtom, 2)
    };
    // The analysed kind's member sites, and which of them are split.
    const memberRefs = [...new Set(payload.templates
        .filter((template) => template.referenceNumbers.length === kind.members && template.edges.length === kind.edges)
        .flatMap((template) => template.referenceNumbers))];
    const split = splitMemberSites(siteRows, memberRefs);
    const splitSites = split.length && 2 * split.length >= memberRefs.length;
    if (split.length) data.split_member_sites = split.map((site) => site.referenceNumber);

    const otherShare = payload.kinds.slice(1).reduce((sum, k) => sum + k.clusters, 0)
        / payload.kinds.reduce((sum, k) => sum + k.clusters, 0);
    if (otherShare > 0.1) {
        checks.push({
            level: 'warn',
            text: `The window gives clusters of different sizes (${payload.kinds.map((k) => `${k.clusters}×${k.members}`).join(', ')}): `
                + 'it may cut through or merge clusters. The statistics are for the most common size only.'
        });
    }
    const ideal = data.average_structure;
    checks.push(payload.ideal.spreadOverSe >= STATIC_SPREAD_SE
        ? {
            level: 'note',
            text: `The average cluster is itself distorted: its edges span ${ideal.edges_A[0]}–${ideal.edges_A.at(-1)} Å, `
                + `${ideal.spread_over_se}× their sampling error (${ideal.edge_se_A} Å) — a static distortion of the average structure.`
        }
        : {
            level: 'note',
            text: `The average cluster is regular within its sampling error: edges ${ideal.edges_A[0]}–${ideal.edges_A.at(-1)} Å `
                + `(± ${ideal.edge_se_A} Å each).`
        });
    // The lattice axes more clusters are elongated or compressed along than
    // by chance: a correlated shape distortion, whatever the edge spread says.
    const favoured = [['long_axis_within_15deg', 'elongated'], ['short_axis_within_15deg', 'compressed']].flatMap(([key, verb]) => (
        Object.entries(data.local[key] || {})
            .filter(([, entry]) => entry.z >= CLUSTER_Z && entry.null > 0 && entry.obs / entry.null >= AXIS_ENHANCEMENT)
            .map(([family, entry]) => ({ family, entry, verb }))
    ));
    const trig = data.local.tetrahedral_modes?.t2_trigonality;
    const trigonalLean = Boolean(trig && significant(trig.z, trig.obs - trig.null, TRIGONALITY_EFFECT));
    const width = data.local.ladder_width_A;
    const widthChanged = significant(width.z, width.change, CLUSTER_EFFECT);
    const breathing = data.local.breathing_sd_A;
    const breathes = significant(breathing.z, breathing.change, CLUSTER_EFFECT);
    const pattern = favoured.length
        ? `more clusters than by chance are ${favoured[0].verb} along ${favoured[0].family} (below)`
        : (trigonalLean ? 'the edge pattern leans to one form (below)' : null);
    const sizeText = '(at each cluster\'s own size)';
    if (kind.edges === 1) {
        checks.push({
            level: 'ok',
            text: `${kind.members === 2 ? 'Dimers have' : 'Each cluster has'} one edge: the shape is the bond length, read from the size `
                + `row (breathing_sd_A, against independent atoms)${pattern ? `; ${pattern}` : ''}.`
        });
    } else if (widthChanged && width.change > 0) {
        checks.push({
            level: 'ok',
            text: `The clusters are distorted beyond independent atoms ${sizeText}: the spread of their edges (longest − shortest, `
                + `${width.obs} Å) is ${Math.round(width.change * 100)} % above the shuffle null (${width.null} Å, z = ${width.z}) — `
                + 'a correlated distortion the data put in.'
        });
    } else if (widthChanged) {
        checks.push({
            level: 'ok',
            text: `The clusters are ${Math.round(-width.change * 100)} % LESS distorted than independent atoms would make them ${sizeText} `
                + `(edge spread ${width.obs} vs ${width.null} Å, z = ${width.z}): the atoms of a cluster move together`
                + (pattern
                    ? `, yet ${pattern}: a correlated shape distortion whose amplitude the independent-atom null already holds (split or disordered sites).`
                    : ', and there is no correlated shape distortion on top.')
        });
    } else if (pattern) {
        checks.push({
            level: 'ok',
            text: `The edge spread (${width.obs} Å) matches the shuffle null ${sizeText} (${width.null} Å; change `
                + `${Math.round((width.change ?? 0) * 100)} %, z = ${width.z}), but the pattern does not: ${pattern}. A correlated shape `
                + 'distortion whose amplitude the sites\' own displacements already hold.'
        });
    } else {
        checks.push({
            level: 'ok',
            text: `No cluster shape distortion beyond independent atoms ${sizeText}: the edge spread (${width.obs} Å) matches the shuffle null `
                + `(${width.null} Å; change ${Math.round((width.change ?? 0) * 100)} %, z = ${width.z}).`
        });
    }
    if (breathes) {
        checks.push({
            level: 'note',
            text: breathing.change > 0
                ? `The clusters breathe: their size (mean edge) fluctuates ${Math.round(breathing.change * 100)} % more than for independent `
                    + `atoms (breathing sd ${breathing.obs} vs ${breathing.null} Å, z = ${breathing.z}) — a collective size fluctuation, `
                    + 'not a shape distortion.'
                : `The cluster size (mean edge) fluctuates ${Math.round(-breathing.change * 100)} % less than for independent atoms `
                    + `(breathing sd ${breathing.obs} vs ${breathing.null} Å, z = ${breathing.z})`
                    + (splitSites
                        ? ': the split sites\' positions combine in a pattern that keeps the size (as a Jahn–Teller distortion does), '
                            + 'not stiff bonds.'
                        : ': the bonds are stiffer than independent atoms make them — the atoms move together.')
        });
    }
    if (trigonalLean) {
        checks.push({
            level: 'note',
            text: trig.obs < trig.null
                ? `The edge pattern leans trigonal: three edges at one vertex differ from the other three (C3v-like; trigonality ${trig.obs} vs ${trig.null}, z = ${trig.z}).`
                : `The edge pattern leans to one edge against its opposite (trigonality ${trig.obs} vs ${trig.null}, z = ${trig.z}).`
        });
    }
    favoured.forEach(({ family, entry, verb }) => {
        checks.push({
            level: 'note',
            text: `More clusters are ${verb} along ${family} than by chance: ${Math.round(entry.obs * 100)} % `
                + `vs ${Math.round(entry.null * 100)} % in the null (z = ${entry.z}).`
        });
    });
    if (splitSites) {
        checks.push({
            level: 'note',
            text: `Split sites: ${split.length} of the ${memberRefs.length} member sites (${split.map((site) => site.referenceNumber).join(', ')}) `
                + 'have two positions (PC1 excess kurtosis ≤ −0.5). The null keeps each site\'s two positions, so it already holds the '
                + 'distortion\'s amplitude: the edge spread cannot show it, and "no distortion" or "less distorted" does not make the '
                + 'clusters regular. Read the pattern (axis families, trigonality), the partial g(r) shell (find_gr_peaks) and the bond angles.'
        });
    }
    if (data.outliers.bonds_outside_window > CLUSTER_OUTLIER_FRACTION) {
        checks.push({
            level: 'warn',
            text: `${Math.round(data.outliers.bonds_outside_window * 100)} % of the cluster bonds lie outside ${lo}–${hi} Å in the configuration: `
                + 'the bond lengths spread wider than the window, or have heavy tails (atoms straying from their cluster).'
        });
    }
    if (data.outliers.contacts_between_clusters_per_atom > CLUSTER_OUTLIER_FRACTION) {
        checks.push({
            level: 'warn',
            text: `${data.outliers.contacts_between_clusters_per_atom} ${element}–${element} contacts per atom inside the window `
                + 'are not cluster bonds: atoms crossing the gap between clusters.'
        });
    }
    // For the answer checks (review L5): what this analysis found.
    const facts = {
        cluster: {
            element,
            members: kind.members,
            window: [lo, hi],
            distorted: (widthChanged && width.change > 0) || Boolean(pattern),
            static: payload.ideal.spreadOverSe >= STATIC_SPREAD_SE,
            breathing: breathes && breathing.change > 0,
            rigid: !splitSites && ((widthChanged && width.change < 0 && !pattern) || (breathes && breathing.change < 0)),
            split: Boolean(splitSites),
            families: favoured.map(({ family, verb }) => `${verb} ${family}`)
        }
    };
    return { data, checks, facts };
};

// --- analyze_modulations ---------------------------------------------------------

// A zone-boundary point is ordered displacement when its power is at least
// this many times the median over all wavevectors (the effect size) AND it
// is significant at SIGNIFICANCE_SIGMA after a family-wise (Šidák)
// correction over every zone-boundary point and element tested. The null:
// for independent displacements a site's normalised power at a half-integer
// q (a real Fourier sum) is Σ_c w_c χ²₁ with Σ w_c = 1, so an element's mean
// over m sites has at least m degrees of freedom: P ≥ χ²_m / m is the
// heaviest tail any site anisotropy allows (conservative). The earlier rule
// (top 1 % of the modes, or rank 1) had no such correction: the validation
// cases found "order" in 6 of 12 configurations without any.
export const MODULATION_RATIO = 3;

// One-sided σ of an upper-tail probability.
const sigmaOf = (p) => (p <= 0 ? Infinity : (p >= 1 ? -Infinity : -normalQuantile(p)));

/**
 * The family-wise significance (one-sided σ) of an element's zone-boundary
 * power `power` (normalised: 1 = independent displacements) over `sites`
 * sites, among `family` points tested.
 */
export const zoneBoundarySigma = (power, sites, family) => {
    if (!(power > 0) || !(sites > 0)) return null;
    const p = regularizedGamma(sites / 2, (sites * power) / 2).upper;
    // Šidák: 1 − (1 − p)^family, computed without cancellation for small p.
    const familyP = -Math.expm1(family * Math.log1p(-Math.min(p, 1 - 1e-16)));
    return sigmaOf(familyP);
};
// A mode this strong, within this fraction of the zone centre, is a
// long-wavelength (box-scale) displacement wave.
export const LONG_WAVE_RATIO = 5;
export const LONG_WAVE_Q = 0.25;

/**
 * The same family-wise significance at a general wavevector (not a
 * zone-boundary point, not a long wave): there the Fourier sum is complex, so
 * for independent displacements a site's normalised power is Σ_c w_c χ²₂/2,
 * and an element's mean over m sites has at least 2m degrees of freedom
 * (P ≥ χ²_2m / 2m); `family` = the compared modes (one per ±q pair). Order
 * there is correlation between cells that does not double the cell:
 * incommensurate, or short-range (the review's L4, approved 2026-10-14).
 */
export const generalModeSigma = (power, sites, family) => {
    if (!(power > 0) || !(sites > 0) || !(family >= 1)) return null;
    const p = regularizedGamma(sites, sites * power).upper;
    const familyP = -Math.expm1(family * Math.log1p(-Math.min(p, 1 - 1e-16)));
    return sigmaOf(familyP);
};
const isZoneBoundary = (q) => q.every((value) => Math.abs(value) < 1e-9 || Math.abs(Math.abs(value) - 0.5) < 1e-9);
const isLongWave = (q) => q.every((value) => Math.abs(value) <= LONG_WAVE_Q + 1e-12);

const qLabel = (q) => `(${q.map((value) => {
    const rounded = Math.round(value * 1000) / 1000;
    return rounded === 0.5 ? '½' : String(rounded + 0);
}).join(', ')})`;

// A mode row is [q label, ratio, percentile, sigma]; sigma reads '>8' past 8σ.
const MODE_ROW = ['string', 'number|null', 'number|null', 'number|string|null'];
const MODULATIONS_SHAPE = {
    note: 'string',
    supercell: ['number', 'number', 'number'],
    wavevectors: 'number',
    elements: [{ element: 'string', sites: 'number', zone_boundary: [MODE_ROW], strongest: [MODE_ROW] }],
    'untestable_zone_boundary?': ['string']
};

export const summarizeModulations = (payload) => {
    const checks = [...parseWarningCheck(payload)];
    const tested = (payload.elements || []).filter((entry) => entry.sites);
    const family = tested.reduce((sum, entry) => sum + entry.zoneBoundary.length, 0);
    const sigmaAt = (entry, point) => (Number.isFinite(point.ratio) && Number.isFinite(entry.medianPower)
        ? zoneBoundarySigma(point.ratio * entry.medianPower, entry.sites, family)
        : null);
    // The general modes: the compared modes over every tested element.
    const generalFamily = tested.reduce((sum, entry) => sum + (Number.isFinite(entry.comparedModes) ? entry.comparedModes : 0), 0);
    const generalSigmaAt = (entry, mode) => (Number.isFinite(mode.ratio) && Number.isFinite(entry.medianPower) && generalFamily > 0
        && !isZoneBoundary(mode.q) && !isLongWave(mode.q)
        ? generalModeSigma(mode.ratio * entry.medianPower, entry.sites, generalFamily)
        : null);
    const data = {
        note: 'zone_boundary rows are [q, ratio, percentile, sigma]: an element\'s displacement power at wavevector q (fractions '
            + 'of the cell\'s reciprocal axes) over its median across all q of the supercell grid (the diffuse background), the '
            + 'share of the compared modes (one per ±q pair, long waves with every |q_i| ≤ 0.25 left out) with less power, and '
            + `the calibrated one-sided significance against independent displacements, corrected over all ${family} zone-boundary `
            + 'points tested (≥ 3 with ratio ≥ 3 is order). strongest rows are [q, ratio, percentile, sigma]: the sigma of a '
            + 'general q (neither zone boundary nor long wave; a complex sum, corrected over the compared modes), null otherwise. '
            + 'A ½ component is order that doubles the cell along that axis; order at a general q is correlation between cells '
            + '(incommensurate or short-range) that does not.',
        supercell: payload.supercell,
        wavevectors: payload.wavevectors,
        elements: tested.map((entry) => ({
            element: entry.element,
            sites: entry.sites,
            // [q, ratio to the median, percentile, σ] — compact rows: 7 points × elements.
            zone_boundary: entry.zoneBoundary.map((point) => {
                const sigma = sigmaAt(entry, point);
                return [qLabel(point.q), round(point.ratio, 2), roundFixed(point.percentile, 1), sigma === Infinity ? '>8' : roundFixed(sigma, 1)];
            }),
            strongest: entry.strongest.slice(0, 3).map((mode) => {
                const sigma = generalSigmaAt(entry, mode);
                return [qLabel(mode.q), round(mode.ratio, 2), roundFixed(mode.percentile, 1), sigma == null ? null : (sigma === Infinity ? '>8' : roundFixed(sigma, 1))];
            })
        }))
    };
    if (payload.untestable?.length) data.untestable_zone_boundary = payload.untestable.map((flags) => qLabel(flags.map((flag) => flag / 2)));
    if (!data.elements.length) {
        checks.push({ level: 'fail', text: 'No site has displacements to analyse (zero spread everywhere).' });
        return { data, checks };
    }
    const ordered = [];
    tested.forEach((entry) => {
        entry.zoneBoundary.forEach((point) => {
            const sigma = sigmaAt(entry, point);
            if (point.ratio >= MODULATION_RATIO && sigma != null && sigma >= SIGNIFICANCE_SIGMA) ordered.push({ element: entry.element, point, sigma });
        });
    });
    const sigmaText = (sigma) => (sigma === Infinity || sigma > 8 ? '> 8σ' : `${round(sigma, 2)}σ`);
    if (ordered.length) {
        ordered.forEach(({ element, point, sigma }) => {
            const axes = ['a', 'b', 'c'].filter((_, axis) => point.q[axis]).join(', ');
            checks.push({
                level: 'ok',
                text: `Ordered ${element} displacements at q = ${qLabel(point.q)}: ${round(point.ratio, 3)}× the median power, `
                    + `${sigmaText(sigma)} against independent displacements (corrected over ${family} zone-boundary points) `
                    + '— a superstructure doubling the cell along '
                    + `${axes}.`
            });
        });
    } else {
        const best = tested
            .flatMap((entry) => entry.zoneBoundary.map((point) => ({ element: entry.element, point, sigma: sigmaAt(entry, point) })))
            .sort((a, b) => (b.sigma ?? -Infinity) - (a.sigma ?? -Infinity))[0];
        checks.push({
            level: 'ok',
            text: 'No order that doubles the cell: no zone-boundary point stands out of the diffuse background'
                + (best ? ` (the strongest, ${best.element} at ${qLabel(best.point.q)}, is ${round(best.point.ratio, 2)}× the median, `
                    + `${best.sigma == null ? 'significance undefined' : sigmaText(best.sigma)} after the correction).` : '.')
        });
    }
    // Correlation between cells at a general wavevector: the strongest modes
    // that are neither zone-boundary points nor long waves, by the same rule.
    const general = [];
    tested.forEach((entry) => {
        entry.strongest.forEach((mode) => {
            const sigma = generalSigmaAt(entry, mode);
            if (sigma != null && mode.ratio >= MODULATION_RATIO && sigma >= SIGNIFICANCE_SIGMA
                && !general.some((item) => item.element === entry.element)) {
                general.push({ element: entry.element, mode, sigma });
            }
        });
    });
    general.forEach(({ element, mode, sigma }) => {
        checks.push({
            level: 'ok',
            text: `Correlated ${element} displacements between cells at q = ${qLabel(mode.q)}: ${round(mode.ratio, 3)}× the median power, `
                + `${sigmaText(sigma)} against independent displacements (corrected over ${generalFamily} compared modes) — order at a `
                + 'general wavevector (incommensurate, or short-range correlation between cells), not a cell doubling.'
        });
    });
    (payload.elements || []).filter((entry) => entry.sites).forEach((entry) => {
        const mode = entry.strongest[0];
        if (mode && mode.ratio >= LONG_WAVE_RATIO && mode.q.every((value) => Math.abs(value) <= LONG_WAVE_Q)) {
            checks.push({
                level: 'note',
                text: `${entry.element}'s strongest mode is a long-wavelength wave, q = ${qLabel(mode.q)} at ${round(mode.ratio, 2)}× the median: `
                    + 'box-scale displacement correlations (common in RMC), not a crystallographic superstructure.'
            });
        }
    });
    if (payload.untestable?.length) {
        checks.push({ level: 'note', text: `An odd supercell axis cannot hold ${data.untestable_zone_boundary.join(', ')}: not tested.` });
    }
    // For the answer check: a superstructure claim needs an ordered point.
    const facts = {
        modulation: {
            elements: tested.map((entry) => entry.element),
            ordered: ordered.map(({ element, point }) => ({ element, q: point.q })),
            general: general.map(({ element, mode }) => ({ element, q: mode.q }))
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
    if (facts) {
        checks.push({
            level: 'note',
            text: `Point group ${describePointGroup(described.pointGroup)}. Take inversion, chirality and polarity from this, not from memory.`
        });
    }
    // For the answer check (llm/copilot/answerAudit.js): the group and what its
    // point group implies.
    const result = { data, checks };
    if (facts && !String(described.spaceGroup).startsWith('≥') && described.spaceGroup !== 'undetermined') {
        result.facts = { groups: [{ symbol: described.spaceGroup, ...facts }] };
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
 * Dashboard's parsed structure; its `basis` (browser-loaded runs only) enables
 * detect_symmetry. Model arguments reach the engines only through the explicit
 * mapping in each run(), never by spreading.
 */
export const buildCopilotTools = ({ request, structure, elements = [] }) => {
    const tools = [];
    const unitCellFrom = (payload) => unitCellVectors(payload?.latticeVectors ?? structure?.latticeVectors, payload?.supercell ?? structure?.supercell);

    if (request) {
        tools.push({
            name: 'list_sites',
            description: 'List every reference site of the RMC configuration with its displacement statistics from a '
                + 'PCA of its displacement cloud: U_iso (Å²), principal RMS amplitudes (Å), anisotropy, Mardia '
                + 'kurtosis, the excess kurtosis along each resolved principal axis, PC1 as a crystal direction [u v w], '
                + 'and for a split site (two positions along PC1) its static displacement with its sampling error. Use it '
                + 'to find the most disordered, anisotropic or anharmonic sites (their `ref`).',
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
                + '(Cartesian and crystal [u v w]) with its calibrated significance, the whole-map test against '
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
                const basis = structure?.basis || [];
                const pooledSites = referenceNumber != null
                    ? 0
                    : (element === 'all' ? basis.length : basis.filter((site) => site.el === element).length);
                return summarizeOrientation(payload, unitCellFrom(payload), { pooledSites });
            }
        });

        tools.push({
            name: 'compute_bond_angles',
            description: 'Bond-angle distribution A–B–C with B the central atom (RMCProfile triplets): every pair of '
                + 'B neighbours, A within [r12_min, r12_max] Å of B and C within [r23_min, r23_max] Å (default: the '
                + 'A–B window). Returns the angle modes, mean ± sd, a 10° histogram, bond lengths and the B '
                + 'coordination, and for about six or four neighbours the split against an ideal octahedron or '
                + 'tetrahedron. Use the whole first coordination shell (find_gr_peaks\' first_shell.window_A): a window '
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
            // not cut a coordination shell before this runs (llm clarify.js).
            shellWindows: (args) => [
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
            ],
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
                return summarizeBondAngles(await request('triplets', params, { signal }), { windowConfirmed: args.window_confirmed === true });
            }
        });

        tools.push({
            name: 'analyze_clusters',
            description: 'Shape of the finite clusters one element forms — the M4 tetrahedra of a lacunar spinel, '
                + 'dimers, trimers — with bonds in [r_min, r_max] Å in the average structure (take the window from '
                + 'find_gr_peaks: first_shell.window_A of that element pair). Every cluster in the configuration is '
                + 'compared with a shuffle null in which each atom keeps its own site\'s displacements but takes them from '
                + 'a random other cluster: is there a distortion beyond independent atoms, which edge pattern '
                + '(tetrahedra: E/T2 modes, trigonal vs one edge), along which lattice axis (<100>/<110>/<111>), and is '
                + 'the average cluster itself distorted. Use it for cluster or polyhedron distortion, Jahn–Teller, '
                + 'dimerization or trimerization questions; single-site displacement directions cannot show these.',
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
                // The site rows tell split member sites; the engine's parse
                // cache makes them cheap, and a failure only drops that note.
                const [payload, sites] = await Promise.all([
                    request('clusters', { element, rMin, rMax }, { signal }),
                    request('sites', {}, { signal }).catch(() => null)
                ]);
                return summarizeClusters(payload, { siteRows: sites?.sites ?? null });
            }
        });

        tools.push({
            name: 'analyze_modulations',
            description: 'Look for order that enlarges the unit cell — a superstructure the folded average structure and '
                + 'detect_symmetry cannot show (cell doubling as in low-temperature phases, octahedral tilts, '
                + 'charge-density waves). Fourier-transforms every site\'s displacements over the supercell\'s cell grid and '
                + 'compares the zone-boundary points (½,0,0)…(½,½,½) and the strongest modes with the configuration\'s own '
                + 'diffuse background (ratio to the median power, percentile). Optionally one element.',
            parameters: {
                type: 'object',
                properties: {
                    element: elementProperty(elements, 'Only this element (default: each element).')
                }
            },
            dataShape: MODULATIONS_SHAPE,
            label: ({ element }) => `Superstructure test${element ? ` · ${element}` : ''}`,
            run: async ({ element }, { signal } = {}) => summarizeModulations(
                await request('modulations', element ? { element } : {}, { signal })
            )
        });
    }

    if (structure?.basis?.length) {
        tools.push({
            name: 'detect_symmetry',
            description: 'Run the space-group finder on the average structure (the folded reference sites) at a '
                + 'position tolerance in Å, as the Detected SG card does: group, number (only when it is certain), '
                + 'operations, and the Wyckoff orbits. Compare tolerances to see which symmetry the local '
                + 'distortion breaks. Does not change the card.',
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

        const floor = symmetryNoiseFloor(structure.basis);
        // The lower group's default tolerance: just above the noise floor.
        const lowerDefault = Number(((floor ?? 0.02) * NOISE_FLOOR_MARGIN).toPrecision(2));
        tools.push({
            name: 'measure_symmetry_lowering',
            description: 'How far each average site sits from the position a higher-symmetry (parent) group puts it at: the '
                + 'size and direction of a coherent shift every copy of a site shares (a polar off-centring of the average '
                + 'structure), which the displacement analyses cannot show — each site\'s cloud is centred on its own mean. The '
                + 'parent is the group at parent_tolerance_A (default 0.2 Å, the Detected SG card\'s); the lower group is the one '
                + `at tolerance_A (default ${lowerDefault} Å, just above the sampling noise floor). Each shift is in Å with its `
                + 'crystal [u v w], measured with the parent on the majority of the atoms, against the noise floor.',
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
    return tools;
};
