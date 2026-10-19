// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// The model-free arms.
//
// Arm A, the default Workbench: what its pages show with their defaults,
// read by fixed rules a competent user would apply (not a straw man): the
// Detected SG card at 0.2 Å, the PCA Ellipsoid table, the Displacement
// Directions page per site (ν = 10, smoothing 2) and the partial g(r) on the
// Dashboard. The cluster and superstructure engines have no page, so arm A
// cannot use them.
//
// Arm A+ (control), the engine sweep: every Copilot engine run exhaustively
// with its own thresholds and fixed rules joining them, no model. What arm B
// adds over A+ is what the language model itself adds. The cluster and
// modulation engines decide for themselves (their verdicts): A+ calls them
// directly and maps their answers to claims. Where the app has a rule, A+
// reads it from where the app decides it (architecture candidate 8,
// SWEEP_RULES below): the split site, the cluster elements and the lower
// group's tolerance.
//
// Each reading is recorded with the value it rests on, so a misleading
// result can be traced to the page that misled.

import { pointGroupFacts } from '../src/llm/context/pointGroups';
import { clusterElementsOf, extractPeaks, shellCatalogue } from '../src/llm/context/shells';
import { lowerGroupTolerance, NOISE_FLOOR_MARGIN, symmetryNoiseFloor } from '../src/llm/context/runContext';
import { describeSymmetry } from '../src/symmetryModel';
import { unitCellVectors } from '../src/pcaCrystalFrame';
import { cartesianToUvw, splitDisplacement } from '../src/copilotTools';
import { isSplitSite, splitSiteEvidence } from '../src/workers/clusters';
import { normalizeVerdict } from './verdict';

// What the arms read of a case: the run as the Copilot is handed it
// (runFolder.js loadRun: a CopilotRun, or the promise loadRun returns), its
// structure, plot files and run settings, its engine, elements and engine
// tools.
const readRun = async (copilotRun) => {
    const run = await copilotRun;
    return {
        structure: run.props.structure,
        plotFiles: run.props.plotFiles,
        runSettings: run.props.runSettings ?? null,
        engine: run.engine,
        elements: run.elements,
        tools: run.tools
    };
};

/** Arm A's reading rules, fixed before any run. */
export const ARM_A_RULES = Object.freeze({
    // PCA Ellipsoid: a split site is a prolate cloud (PC1 rms over PC2 rms)
    // whose PC1 excess kurtosis is clearly negative (two maxima); a Gaussian
    // ellipsoid, however anisotropic, is thermal motion.
    PROLATE_RATIO: 1.3,
    SPLIT_KURTOSIS: -0.5,
    // Displacement Directions page defaults, and its significance bar.
    FREQUENCY: 10,
    SMOOTHING: 2,
    SIGMA: 3,
    // Partial g(r): maxima below this (the bonded range, as the eye reads
    // the plot) at least this far apart are a split shell, unless the
    // average structure itself has distances that far apart there.
    BOND_MAX_A: 3.2,
    SPLIT_SHELL_A: 0.15,
    // Detected SG card default.
    CARD_TOLERANCE_A: 0.2
});

const round = (value, digits = 3) => (Number.isFinite(value) ? Number(value.toPrecision(digits)) : null);

// The crystal [u v w] of a Cartesian direction, cell from the run.
const uvwOf = (structure, vector) => cartesianToUvw(vector, unitCellVectors(structure.latticeVectors, structure.supercell));

const symmetryReading = (described, tolerance, source) => {
    const facts = described?.pointGroup ? pointGroupFacts(described.pointGroup) : null;
    const certain = described && described.spaceGroup !== 'undetermined' && !String(described.spaceGroup).startsWith('≥');
    return {
        reading: { page: source, what: `space group at ${tolerance} Å`, value: described?.spaceGroup ?? null, point_group: described?.pointGroup ?? null },
        verdict: {
            space_group: certain ? described.spaceGroup : null,
            chiral: facts ? facts.sohncke : null,
            centrosymmetric: facts ? facts.centrosymmetric : null,
            polar: facts ? facts.polar : null
        }
    };
};

// PCA Ellipsoid rows → split sites (one claim per element: its largest).
const splitSites = async (run, readings) => {
    const payload = await run.engine('sites');
    const found = new Map();
    for (const site of payload.sites || []) {
        if (!site.rms || !site.axes) continue;
        const ratio = site.rms[0] / site.rms[1];
        const kurtosis = site.excessKurtosis?.[0];
        const split = ratio >= ARM_A_RULES.PROLATE_RATIO && site.axisResolved?.[0] && Number.isFinite(kurtosis) && kurtosis <= ARM_A_RULES.SPLIT_KURTOSIS;
        const magnitude = Math.sqrt(Math.max(site.rms[0] ** 2 - site.rms[1] ** 2, 0));
        readings.push({
            page: 'PCA Ellipsoid', what: `site ${site.referenceNumber} ${site.element}`,
            value: { rms_A: site.rms.map((r) => round(r)), pc1_over_pc2: round(ratio), pc1_excess_kurtosis: round(kurtosis), pc1_uvw: uvwOf(run.structure, site.axes[0]) },
            reading: split ? `split site along PC1, static shift ≈ ${round(magnitude)} Å` : 'thermal ellipsoid'
        });
        if (!split) continue;
        const current = found.get(site.element);
        if (!current || magnitude > current.magnitude) found.set(site.element, { site, magnitude });
    }
    return { payload, found };
};

// Displacement Directions page on one site, at the page's defaults.
const directionsOf = async (run, site, readings) => {
    const payload = await run.engine('orientation', {
        referenceNumber: site.referenceNumber, frequency: ARM_A_RULES.FREQUENCY, smoothing: ARM_A_RULES.SMOOTHING, geometry: false
    });
    const significant = payload.peakSignificance >= ARM_A_RULES.SIGMA;
    const uvw = Array.isArray(payload.peakDirection) ? uvwOf(run.structure, payload.peakDirection) : null;
    readings.push({
        page: 'Displacement Directions', what: `site ${site.referenceNumber} ${site.element}`,
        value: { peak_sigma: round(payload.peakSignificance), enhancement: round(payload.peakEnhancement), peak_uvw: uvw, map_sigma: round(payload.mapSignificance) },
        reading: significant ? `preferred direction [${uvw?.join(' ')}]` : 'no significant direction'
    });
    return significant ? uvw : null;
};

// The average structure's distinct A–B distances up to rMax (Å).
const averageDistances = (structure, a, b, rMax) => {
    const cell = unitCellVectors(structure.latticeVectors, structure.supercell);
    const cart = (f) => [0, 1, 2].map((k) => f[0] * cell[0][k] + f[1] * cell[1][k] + f[2] * cell[2][k]);
    const out = [];
    const from = structure.basis.filter((site) => site.el === a);
    const to = structure.basis.filter((site) => site.el === b);
    for (const s of from) {
        for (const t of to) {
            for (let i = -1; i <= 1; i += 1) {
                for (let j = -1; j <= 1; j += 1) {
                    for (let k = -1; k <= 1; k += 1) {
                        const d = cart([t.frac[0] + i - s.frac[0], t.frac[1] + j - s.frac[1], t.frac[2] + k - s.frac[2]]);
                        const r = Math.hypot(...d);
                        if (r > 1e-6 && r <= rMax) out.push(r);
                    }
                }
            }
        }
    }
    return out.sort((x, y) => x - y);
};

// Dashboard partial g(r): split bonded shells the average structure does not explain.
const splitShells = (run, readings) => {
    const partials = run.plotFiles.find((file) => file.plotKind === 'pdf_partials');
    const found = [];
    for (const series of partials?.plotData?.series || []) {
        const match = /^([A-Z][a-z]?)\s*-\s*([A-Z][a-z]?)$/.exec((series.label || '').trim());
        if (!match || match[1] === match[2]) continue;
        const peaks = extractPeaks(series.x, series.y, { rMax: ARM_A_RULES.BOND_MAX_A, maxPeaks: 4 }).map((peak) => peak.r);
        const width = peaks.length > 1 ? peaks.at(-1) - peaks[0] : 0;
        if (width < ARM_A_RULES.SPLIT_SHELL_A) continue;
        // The less abundant element is the polyhedron's centre.
        const [a, b] = [match[1], match[2]];
        const centre = (run.structure.elementCounts[a] || 0) <= (run.structure.elementCounts[b] || 0) ? a : b;
        const partner = centre === a ? b : a;
        const average = averageDistances(run.structure, centre, partner, peaks.at(-1) + 0.3).filter((r) => r >= peaks[0] - 0.3);
        const averageWidth = average.length ? average.at(-1) - average[0] : 0;
        const explained = averageWidth >= ARM_A_RULES.SPLIT_SHELL_A;
        readings.push({
            page: 'Dashboard partial g(r)', what: `${series.label} first shell`,
            value: { peaks_A: peaks.map((r) => round(r)), average_structure_spread_A: round(averageWidth) },
            reading: explained ? 'split, as in the average structure' : `split shell around ${centre}: ${round(width)} Å between its maxima`
        });
        if (!explained) found.push({ centre, partner, width });
    }
    return found;
};

/** Arm A: the default pages, read by ARM_A_RULES. `copilotRun`: loadRun's run (or its promise). */
export const runArmA = async (copilotRun) => {
    const run = await readRun(copilotRun);
    const readings = [];
    const card = describeSymmetry(run.structure, ARM_A_RULES.CARD_TOLERANCE_A);
    const symmetry = symmetryReading(card, ARM_A_RULES.CARD_TOLERANCE_A, 'Detected SG card');
    readings.push(symmetry.reading);
    const { found } = await splitSites(run, readings);
    const features = [];
    for (const [element, { site, magnitude }] of found) {
        const peak = await directionsOf(run, site, readings);
        features.push({ kind: 'off_centre', element, direction_uvw: peak ?? uvwOf(run.structure, site.axes[0]), magnitude_A: round(magnitude) });
    }
    for (const shell of splitShells(run, readings)) {
        // An off-centred centre splits its own shell: that explains it.
        if (found.has(shell.centre)) continue;
        features.push({ kind: 'jahn_teller', element: shell.centre, magnitude_A: round(shell.width) });
    }
    return { arm: 'A', readings, verdict: normalizeVerdict({ features, ...symmetry.verdict }) };
};

/**
 * Arm A+'s rules (the engine sweep). Since architecture candidate 8
 * (2026-10-19) every one is the app's, imported from the module that decides
 * it, so A+ and the Copilot cannot read a run two ways:
 *
 *   split site      isSplitSite (workers/clusters.js: PC1 resolved, κ₁ ≤ −0.5
 *                   and ≥ 3 sampling errors below 0), with list_sites' static
 *                   displacement δ = √(rms1² − rms⊥²) (splitDisplacement);
 *                   until then A+ used arm A's (PC1/PC2 ≥ 1.3, κ₁ ≤ −0.5).
 *   cluster elements  clusterElementsOf (llm/context/shells.js): an element
 *                   whose own first shell holds 0.5–4.5 neighbours, the
 *                   evidence gate's rule; until then a same-element first
 *                   shell peaking below 3.6 Å.
 *   lower group     lowerGroupTolerance (llm/context/runContext.js):
 *                   NOISE_FLOOR_MARGIN (1.25) × the sampling noise floor,
 *                   measure_symmetry_lowering's default and the Workbench's
 *                   tolerance question's margin; until then the floor itself.
 *
 *   polarisation    the modulation engine's (plan item 4.3, protocol
 *                   2026-10-19.12): the leading axis of the element's power
 *                   tensor at the q that decided the order, its [u v w];
 *                   until then the PC1 of the element's most prolate site
 *                   (PC1/PC2 ≥ ARM_A_RULES.PROLATE_RATIO), a rule of A+'s own.
 *
 * Arm A's rules (ARM_A_RULES) are the frozen protocol's and stay as they were.
 * What remains A+'s own is only how it joins the engines: one claim per
 * element, an ordered mode explains its element's split sites, and an
 * off-centred centre its own split shell.
 *
 * The table holds the app's functions themselves, and A+ calls them through
 * it: a test checks they are the app's (one definition).
 */
export const SWEEP_RULES = Object.freeze({
    splitSite: isSplitSite,
    staticDisplacement: splitDisplacement,
    clusterElements: clusterElementsOf,
    lowerGroupTolerance
});

/**
 * A+'s split sites, by the app's rule: each site row of the 'sites' request
 * read as list_sites reads it (isSplitSite, through its evidence record: κ₁,
 * its z = κ₁/√(24/n), PC1 resolved), a split site's static displacement
 * δ ± its sampling error as list_sites gives it (delta_A, delta_se_A), and
 * per element the split site with the largest δ (one claim per element).
 * `run` holds the engine and the structure (the cell for [u v w]).
 */
export const appSplitSites = async (run, readings) => {
    const payload = await run.engine('sites');
    const found = new Map();
    for (const site of payload.sites || []) {
        if (!site.rms || !site.axes) continue;
        // The values the rule reads, for the record; the rule decides.
        const evidence = splitSiteEvidence(site);
        const split = SWEEP_RULES.splitSite(site) ? SWEEP_RULES.staticDisplacement(site) : null;
        readings.push({
            page: 'list_sites', what: `site ${site.referenceNumber} ${site.element}`,
            value: {
                rms_A: site.rms.map((r) => round(r)),
                pc1_resolved: Boolean(site.axisResolved?.[0]),
                pc1_excess_kurtosis: round(evidence.observed),
                pc1_kurtosis_z: round(evidence.z),
                pc1_uvw: uvwOf(run.structure, site.axes[0]),
                ...(split ? { delta_A: round(split.delta), delta_se_A: round(split.se, 2) } : {})
            },
            reading: split ? `split site along PC1, static shift δ = ${round(split.delta)} ± ${round(split.se, 2)} Å` : 'not a split site'
        });
        if (!split) continue;
        const current = found.get(site.element);
        if (!current || split.delta > current.magnitude) found.set(site.element, { site, magnitude: split.delta });
    }
    return { payload, found };
};

// An engine Finding as the readings record it: "answer (because)".
const said = (finding) => `${finding.answer} (${finding.because})`;
const qText = (q) => `(${q.map((value) => round(value, 3) ?? 0).join(', ')})`;
const sameQ = (a, b) => a.length === b.length && a.every((value, axis) => value === b[axis]);
// A direct-lattice [u v w] as A+'s other direction claims give it (two decimals, cartesianToUvw's).
const uvwText = (uvw) => uvw.map((value) => Number(value.toFixed(2)) + 0);

// The point of an element's payload entry at q (a zone-boundary point or a strongest mode).
const pointAt = (entry, q) => [...(entry.zoneBoundary || []), ...(entry.strongest || [])].find((point) => sameQ(point.q, q)) ?? null;

// The modulation engine's verdicts A+ reads, per element tested, as a
// reading, with the polarisation of each ordered point (the direction a
// claim rests on). Diffuse order between cells (verdict.diffuseOrder, plan
// item 4.4) is not among them: the verdict has no feature kind for it, and
// adding one changes the scoring (the user's call).
const MODULATION_FINDINGS_READ = ['cellDoubling', 'generalOrder', 'longWave'];
const modulationReading = (payload) => {
    const tested = (payload.elements || []).filter((entry) => entry.sites);
    const named = (name) => tested.filter((entry) => entry.verdict[name].answer === true)
        .map((entry) => `${entry.element} at ${entry.verdict[name].decidedBy.map(qText).join(', ')}`);
    const polarisations = (entry) => ['cellDoubling', 'generalOrder']
        .flatMap((name) => (entry.verdict[name].answer === true ? entry.verdict[name].decidedBy : []))
        .map((q) => ({ q, polarisation: pointAt(entry, q)?.polarisation }))
        .filter((item) => item.polarisation)
        .map(({ q, polarisation }) => ({ q: qText(q), uvw: uvwText(polarisation.uvw), share: round(polarisation.share) }));
    return {
        value: Object.fromEntries(tested.map((entry) => {
            const verdicts = Object.fromEntries(MODULATION_FINDINGS_READ.map((name) => [name, said(entry.verdict[name])]));
            const rows = polarisations(entry);
            return [entry.element, rows.length ? { ...verdicts, polarisation: rows } : verdicts];
        })),
        reading: [['cellDoubling', 'cell doubling'], ['generalOrder', 'order at a general q'], ['longWave', 'long wave']]
            .map(([name, label]) => `${label}: ${named(name).join('; ') || 'none'}`).join('; ')
    };
};

// The cluster engine's verdict as a reading: its answers, and the two rows they mostly rest on.
const clusterReading = (payload) => {
    const { verdict } = payload;
    const row = (record) => (record ? { obs: round(record.observed), null: round(record.null), z: round(record.z, 2), change: round(record.change, 2) } : null);
    if (verdict.window.extendedNetwork.because !== 'finite_clusters') {
        return { value: { window: said(verdict.window.extendedNetwork) }, reading: `no clusters: ${verdict.window.extendedNetwork.because}` };
    }
    const names = ['distorted', 'lessDistorted', 'breathing', 'rigid', 'staticDistortion', 'splitMemberSites', 'trigonalLean', 'favouredAxes'];
    return {
        value: {
            ...Object.fromEntries(names.map((name) => [name, said(verdict[name])])),
            ladder_width: row(verdict.lessDistorted.evidence[0]),
            breathing_sd: row(verdict.breathing.evidence[0])
        },
        reading: names.map((name) => `${name} ${said(verdict[name])}`).join('; ')
    };
};

/** Arm A+: every engine, its own thresholds, fixed joining rules; no model. `copilotRun` as in runArmA. */
export const runArmSweep = async (copilotRun) => {
    const run = await readRun(copilotRun);
    const readings = [];
    // The run's engine tools, as the Copilot is offered them.
    const tools = Object.fromEntries(run.tools.map((tool) => [tool.name, tool]));
    // Symmetry at the lower group's tolerance, the app's: NOISE_FLOOR_MARGIN ×
    // the sampling noise floor (a group change below the floor is noise, and
    // a group counts only where it holds clearly above it).
    const floor = symmetryNoiseFloor(run.structure.basis);
    const tolerance = SWEEP_RULES.lowerGroupTolerance(run.structure.basis);
    const symmetry = symmetryReading(describeSymmetry(run.structure, tolerance), tolerance, `detect_symmetry at ${NOISE_FLOOR_MARGIN} × the noise floor`);
    symmetry.reading.noise_floor_A = round(floor);
    readings.push(symmetry.reading);

    const { found } = await appSplitSites(run, readings);
    const features = [];

    // Superstructure: the modulation engine's verdict, a cell doubling (an
    // ordered zone-boundary point) or (batch 1, review L4) order at a general
    // wavevector, each a modulation claim at the q that decided it.
    const modulations = await run.engine('modulations', {});
    readings.push({ page: 'analyze_modulations', what: 'zone-boundary points and general modes (engine verdict)', ...modulationReading(modulations) });
    const tested = (modulations.elements || []).filter((entry) => entry.sites);
    const orderedAt = (name) => tested.filter((entry) => entry.verdict[name].answer === true)
        .map((entry) => ({ entry, q: entry.verdict[name].decidedBy[0] }));
    const modulated = new Set();
    for (const { entry, q } of [...orderedAt('cellDoubling'), ...orderedAt('generalOrder')]) {
        if (modulated.has(entry.element)) continue;
        modulated.add(entry.element);
        // Polarisation: the engine's, at the q that decided the order (plan item 4.3).
        const polarisation = pointAt(entry, q)?.polarisation;
        features.push({ kind: 'modulation', element: entry.element, q, direction_uvw: polarisation ? uvwText(polarisation.uvw) : null });
    }

    // Split sites not explained by an ordered mode.
    for (const [element, { site, magnitude }] of found) {
        if (modulated.has(element)) continue;
        const peak = await directionsOf(run, site, readings);
        features.push({ kind: 'off_centre', element, direction_uvw: peak ?? uvwOf(run.structure, site.axes[0]), magnitude_A: round(magnitude) });
    }
    for (const shell of splitShells(run, readings)) {
        if (found.has(shell.centre)) continue;
        features.push({ kind: 'jahn_teller', element: shell.centre, magnitude_A: round(shell.width) });
    }

    // A lowering of the average structure (batch 1, review L4): when the
    // group at the lower tolerance is not the card's, the shift tool measures
    // each element's coherent shift from the card's group, above the floor
    // (its lower group at the same tolerance, the tool's default).
    if (tools.measure_symmetry_lowering && symmetry.verdict.space_group
        && symmetry.verdict.space_group !== describeSymmetry(run.structure, ARM_A_RULES.CARD_TOLERANCE_A)?.spaceGroup) {
        const lowering = await tools.measure_symmetry_lowering.run({ parent_tolerance_A: ARM_A_RULES.CARD_TOLERANCE_A, tolerance_A: tolerance });
        readings.push({ page: 'measure_symmetry_lowering', what: `shifts from the ${ARM_A_RULES.CARD_TOLERANCE_A} Å group`, value: lowering.facts?.lowering ?? null, reading: lowering.checks.map((c) => c.text).join(' ') });
        for (const shifted of lowering.facts?.lowering?.shifted || []) {
            features.push({ kind: 'symmetry_lowering', element: shifted.element, direction_uvw: shifted.direction_uvw, magnitude_A: shifted.shift_A });
        }
    }

    // Clusters: every element whose own first shell holds 0.5–4.5 neighbours
    // (the evidence gate's rule, clusterElementsOf), read by the cluster
    // engine's verdict: distorted (the edge spread above the shape null, or a
    // pattern: an axis family, a trigonal lean), breathing, rigid (never with
    // split member sites, where the null mixes the sites' static choices). The
    // shell and its window are the run's one shell reading (shells.js), the
    // window find_gr_peaks gives the model.
    const CLUSTER_CLAIMS = [['distorted', 'cluster_distortion'], ['breathing', 'cluster_breathing'], ['rigid', 'rigid_cluster']];
    const catalogue = shellCatalogue({ plotFiles: run.plotFiles, structure: run.structure, runSettings: run.runSettings });
    for (const element of SWEEP_RULES.clusterElements(catalogue)) {
        const { window } = catalogue.shell(element, element);
        const result = await run.engine('clusters', { element, rMin: window[0], rMax: window[1] });
        readings.push({ page: 'analyze_clusters', what: `${element} clusters in ${window[0]}–${window[1]} Å (engine verdict)`, ...clusterReading(result) });
        CLUSTER_CLAIMS.forEach(([finding, kind]) => {
            if (result.verdict[finding].answer === true) features.push({ kind, element });
        });
    }
    return { arm: 'A+', readings, verdict: normalizeVerdict({ features, ...symmetry.verdict }) };
};
