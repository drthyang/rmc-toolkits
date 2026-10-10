// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// The AI Copilot's analysis tools: the host side of the llm module's tool
// contract (src/llm/copilot/tools.js). Each tool runs one of the Workbench
// engines the pages use — the shared PCA worker for a browser-loaded run, the
// Flask /api/pca/* and /api/triplets routes for a server folder (both through
// useSiteCloud's requestAnalysis), and the symmetry finder — and returns a
// compact summary for the model plus deterministic checks on the full result.
// The tools only compute: none changes a page, the shared tolerance, or a file.
// The summarizers are pure and exported for tests.

import { useEffect, useMemo, useRef } from 'react';
import { isStaticMode } from './browserData';
import { requestAnalysis } from './useSiteCloud';
import { describeSymmetry, orbitLabel } from './symmetryModel';
import { siteLabel } from './siteLabel';
import { crystalOrientationRows, invert3, matVec3, transpose3, unitCellVectors } from './pcaCrystalFrame';

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

// --- list_sites ------------------------------------------------------------------

// Rows are ranked by |kurtosis|, then U_iso — the order the run context uses —
// so trimming to the result budget drops the least remarkable sites first.
export const summarizeSites = (payload, { element = null } = {}) => {
    const all = payload?.sites || [];
    const picked = element ? all.filter((site) => site.element === element || site.elementCounts?.[element]) : all;
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
        if (site.zeroSpread) row.zero_spread = true;
        return row;
    }).sort((a, b) => (
        Math.abs(b.kurtosis ?? 0) - Math.abs(a.kurtosis ?? 0) || (b.U_iso_A2 ?? -Infinity) - (a.U_iso_A2 ?? -Infinity)
    ));
    const ranked = (key) => [...sites]
        .filter((site) => Number.isFinite(site[key]))
        .sort((a, b) => Math.abs(b[key]) - Math.abs(a[key]))
        .slice(0, 3)
        .map((site) => site.ref);
    const data = {
        note: 'kurtosis = Mardia multivariate excess kurtosis (0 Gaussian; >0 peaked/heavy-tailed; <0 flat or split); '
            + 'anisotropy = rms1/rms3; sites sorted by |kurtosis|, largest first',
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

export const summarizeOrientation = (payload, unitCell) => {
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
        checks.push({ level: 'ok', text: `The map departs from isotropy at ${round(payload.mapSignificance)}σ.` });
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
    return { data, checks };
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

// The angle distribution as fractions in 18 bins of 10°.
const coarseHistogram = (centers, counts, total) => {
    const bins = new Array(18).fill(0);
    centers.forEach((center, index) => {
        bins[Math.min(17, Math.floor(center / 10))] += counts[index];
    });
    return bins.map((count) => (total ? round(count / total, 2) : 0));
};

const lengthSummary = (lengths) => (lengths ? {
    mean_A: round(lengths.meanLength, 4),
    unique_bonds: lengths.uniqueBonds
} : null);

export const summarizeBondAngles = (payload) => {
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
    } else {
        const unbonded = apexTotal ? (coordination[0] || 0) / apexTotal : 0;
        if (unbonded > UNBONDED_APEX_FRACTION) {
            checks.push({
                level: 'warn',
                text: `${Math.round(unbonded * 100)} % of the ${b} atoms have no ${a}–${b} bond in the window: it misses most of the shell.`
            });
        } else if (apexTotal) {
            const modalCount = Math.max(...coordination);
            const modal = coordination.indexOf(modalCount);
            const share = modalCount / apexTotal;
            const present = coordination.map((count, bonds) => (count > 0 ? bonds : null)).filter((bonds) => bonds != null);
            if (share < MODAL_COORDINATION_SHARE) {
                checks.push({
                    level: 'warn',
                    text: `Only ${Math.round(share * 100)} % of the ${b} atoms have the most common count (${modal} ${a}–${b} bonds); `
                        + `the counts spread from ${present[0]} to ${present.at(-1)}. `
                        + 'The window may cut through a coordination shell: compare it with the g(r) peaks (find_gr_peaks) '
                        + 'and widen it to the shell\'s outer edge before reading the angles or the coordination.'
                });
            }
        }
        if (!checks.some((check) => check.level === 'warn')) {
            checks.push({ level: 'ok', text: `${total} angles around ${apexTotal} ${b} atoms (mean coordination ${round(meanCoordination)}).` });
        }
    }
    return { data, checks };
};

// --- detect_symmetry -------------------------------------------------------------

export const summarizeSymmetry = (described, tolerance) => {
    if (!described) {
        return {
            data: { tolerance_A: tolerance },
            checks: [{ level: 'fail', text: 'The structure has no average sites to analyse.' }]
        };
    }
    const data = {
        tolerance_A: tolerance,
        space_group: described.spaceGroup,
        number: described.spaceGroupNumber ?? null,
        point_group: described.pointGroup,
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
    return { data, checks };
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
                + 'PCA of its displacement cloud: U_iso (Å²), principal RMS amplitudes (Å), anisotropy and Mardia '
                + 'kurtosis. Use it to find the most disordered, anisotropic or anharmonic sites (their `ref`).',
            parameters: {
                type: 'object',
                properties: {
                    element: elementProperty(elements, 'Only the sites holding this element.')
                }
            },
            label: ({ element }) => `Site displacement table${element ? ` · ${element}` : ''}`,
            page: 'ellipsoids',
            run: async ({ element }, { signal } = {}) => summarizeSites(await request('sites', {}, { signal }), { element })
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
            label: ({ site_ref: ref, element }) => `Displacement directions ${ref != null ? `site ${ref}` : (element || '')}`.trim(),
            page: 'orientation',
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
                return summarizeOrientation(payload, unitCellFrom(payload));
            }
        });

        tools.push({
            name: 'compute_bond_angles',
            description: 'Bond-angle distribution A–B–C with B the central atom (RMCProfile triplets): every pair of '
                + 'B neighbours, A within [r12_min, r12_max] Å of B and C within [r23_min, r23_max] Å (default: the '
                + 'A–B window). Returns the angle modes, mean ± sd, a 10° histogram, bond lengths and the B '
                + 'coordination. Pick windows around the first g(r) peak (find_gr_peaks).',
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
                    bin_width: { type: 'number', minimum: 0.05, maximum: 30, description: 'Histogram bin (degrees). Default 1.' }
                },
                required: ['end1', 'apex', 'end2', 'r12_min', 'r12_max']
            },
            label: ({ end1, apex, end2, r12_min: lo, r12_max: hi }) => (
                `Bond angles ${end1}–${apex}–${end2}${Number.isFinite(lo) && Number.isFinite(hi) ? ` · ${lo}–${hi} Å` : ''}`
            ),
            page: 'geometry',
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
                return summarizeBondAngles(await request('triplets', params, { signal }));
            }
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
            label: ({ tolerance_A: tol }) => `Symmetry at ${tol} Å`,
            page: 'dashboard',
            run: async ({ tolerance_A: tol }) => summarizeSymmetry(describeSymmetry(structure, tol), tol)
        });
    }
    return tools;
};

/**
 * The Copilot's tools for the open run. A browser-loaded run (a picked folder
 * or the Demo) is analysed from its `.rmc6f` text, read lazily on the first
 * request and cached per file, in a worker of the Copilot's own: a stopped or
 * timed-out request terminates it (the next request starts a fresh one), and
 * the pages' shared worker never queues behind the Copilot. A typed server
 * folder goes to Flask, the request cancelled on abort. The list is rebuilt
 * only when the run, the folder or the structure changes.
 */
export default function useCopilotTools({ directory, localRun, structure }) {
    const textRef = useRef({ file: null, text: null });
    const workerRef = useRef(null);
    const localFile = localRun?.structureFile?.sourceFile || null;
    const canRequest = Boolean(localFile) || (!localRun && !isStaticMode());

    useEffect(() => () => {
        workerRef.current?.terminate();
        workerRef.current = null;
    }, []);

    return useMemo(() => {
        const readText = async () => {
            if (textRef.current.file !== localFile) {
                textRef.current = { file: localFile, text: await localFile.text() };
            }
            return textRef.current.text;
        };
        const workerRequest = async (kind, params, signal) => {
            const text = await readText();
            if (!workerRef.current) {
                workerRef.current = new Worker(new URL('./workers/pcaKdeWorker.js', import.meta.url), { type: 'module' });
            }
            const worker = workerRef.current;
            const onAbort = () => {
                if (workerRef.current === worker) workerRef.current = null;
                worker.terminate();
            };
            signal?.addEventListener('abort', onAbort, { once: true });
            try {
                return await requestAnalysis(kind, params, { rmc6fText: text, worker, signal });
            } finally {
                signal?.removeEventListener('abort', onAbort);
            }
        };
        const request = canRequest
            ? (kind, params, { signal } = {}) => (localFile
                ? workerRequest(kind, params, signal)
                : requestAnalysis(kind, params, { directory, signal }))
            : null;
        const elements = Object.keys(structure?.elementCounts || {}).sort();
        return buildCopilotTools({ request, structure, elements });
    }, [canRequest, localFile, directory, structure]);
}
