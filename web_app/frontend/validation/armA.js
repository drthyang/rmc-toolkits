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
// adds over A+ is what the language model itself adds.
//
// Each reading is recorded with the value it rests on, so a misleading
// result can be traced to the page that misled.

import { pointGroupFacts } from '../src/llm/context/pointGroups';
import { extractPeaks, firstShell } from '../src/llm/context/pairCorrelations';
import { symmetryNoiseFloor } from '../src/llm/context/runContext';
import { describeSymmetry } from '../src/symmetryModel';
import { unitCellVectors } from '../src/pcaCrystalFrame';
import {
    AXIS_ENHANCEMENT, CLUSTER_EFFECT, CLUSTER_Z, cartesianToUvw
} from '../src/copilotTools';
import { normalizeVerdict } from './verdict';

// What the arms read of a case: the run as the Copilot is handed it
// (runFolder.js loadRun: a CopilotRun, or the promise loadRun returns), its
// structure and plot files, its engine, elements and engine tools.
const readRun = async (copilotRun) => {
    const run = await copilotRun;
    return { structure: run.props.structure, plotFiles: run.props.plotFiles, engine: run.engine, elements: run.elements, tools: run.tools };
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

/** Arm A+'s extra rules (the engine sweep), fixed before any run. */
export const SWEEP_RULES = Object.freeze({
    // Same-element first shells shorter than this are tried as clusters.
    CLUSTER_SHELL_MAX_A: 3.6
});

const significantChange = (row) => Number.isFinite(row?.z) && Math.abs(row.z) >= CLUSTER_Z
    && Number.isFinite(row?.change) && Math.abs(row.change) >= CLUSTER_EFFECT;

/** Arm A+: every engine, its own thresholds, fixed joining rules; no model. `copilotRun` as in runArmA. */
export const runArmSweep = async (copilotRun) => {
    const run = await readRun(copilotRun);
    const readings = [];
    // The run's engine tools, as the Copilot is offered them.
    const tools = Object.fromEntries(run.tools.map((tool) => [tool.name, tool]));
    // Symmetry at the sampling noise floor: the local-symmetry skill's rule
    // (a group change below the floor is noise).
    const floor = symmetryNoiseFloor(run.structure.basis);
    const tolerance = Number((floor ?? ARM_A_RULES.CARD_TOLERANCE_A).toPrecision(2));
    const symmetry = symmetryReading(describeSymmetry(run.structure, tolerance), tolerance, 'detect_symmetry at the noise floor');
    readings.push(symmetry.reading);

    const { payload: sites, found } = await splitSites(run, readings);
    const features = [];

    // Superstructure: the modulation tool's own rules, for an ordered
    // zone-boundary point and (batch 1, review L4) a general wavevector.
    const modulations = await tools.analyze_modulations.run({});
    readings.push({ page: 'analyze_modulations', what: 'zone-boundary points and general modes', value: null, reading: modulations.checks.map((c) => c.text).join(' ') });
    const modulated = new Set();
    for (const { element, q } of [...modulations.facts.modulation.ordered, ...(modulations.facts.modulation.general || [])]) {
        if (modulated.has(element)) continue;
        modulated.add(element);
        // Polarisation: PC1 of the element's most prolate site.
        const site = (sites.sites || []).filter((s) => s.element === element && s.axes)
            .sort((x, y) => (y.rms[0] / y.rms[1]) - (x.rms[0] / x.rms[1]))[0];
        const prolate = site && site.rms[0] / site.rms[1] >= ARM_A_RULES.PROLATE_RATIO;
        features.push({ kind: 'modulation', element, q, direction_uvw: prolate ? uvwOf(run.structure, site.axes[0]) : null });
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
    // group at the noise floor is not the card's, the shift tool measures
    // each element's coherent shift from the card's group, above the floor.
    if (tools.measure_symmetry_lowering && symmetry.verdict.space_group
        && symmetry.verdict.space_group !== describeSymmetry(run.structure, ARM_A_RULES.CARD_TOLERANCE_A)?.spaceGroup) {
        const lowering = await tools.measure_symmetry_lowering.run({ parent_tolerance_A: ARM_A_RULES.CARD_TOLERANCE_A, tolerance_A: tolerance });
        readings.push({ page: 'measure_symmetry_lowering', what: `shifts from the ${ARM_A_RULES.CARD_TOLERANCE_A} Å group`, value: lowering.facts?.lowering ?? null, reading: lowering.checks.map((c) => c.text).join(' ') });
        for (const shifted of lowering.facts?.lowering?.shifted || []) {
            features.push({ kind: 'symmetry_lowering', element: shifted.element, direction_uvw: shifted.direction_uvw, magnitude_A: shifted.shift_A });
        }
    }

    // Clusters: every element whose own first shell is short and closed.
    const partials = run.plotFiles.find((file) => file.plotKind === 'pdf_partials');
    for (const series of partials?.plotData?.series || []) {
        const match = /^([A-Z][a-z]?)\s*-\s*([A-Z][a-z]?)$/.exec((series.label || '').trim());
        if (!match || match[1] !== match[2]) continue;
        const shell = firstShell(series.x, series.y);
        if (!shell || shell.boundaryR == null || shell.peakR > SWEEP_RULES.CLUSTER_SHELL_MAX_A) continue;
        const element = match[1];
        const window = [Number(shell.onsetR.toFixed(2)), Number(shell.boundaryR.toFixed(2))];
        const result = await tools.analyze_clusters.run({ element, r_min: window[0], r_max: window[1] });
        const local = result.data.local;
        readings.push({
            page: 'analyze_clusters', what: `${element} clusters in ${window[0]}–${window[1]} Å`,
            value: local ? { ladder_width: local.ladder_width_A, breathing_sd: local.breathing_sd_A } : null,
            reading: result.checks.map((c) => c.text).join(' ')
        });
        if (!local) continue;
        const width = local.ladder_width_A;
        const breathing = local.breathing_sd_A;
        const kinds = new Set();
        const breathes = significantChange(breathing) && breathing.change > 0;
        // "Less than independent atoms" means rigid only for unimodal sites:
        // for split sites the null mixes static choices, and a collective
        // breathing alone lowers the shape spread.
        const rigidReadable = !found.has(element) && !breathes;
        if (significantChange(width)) {
            if (width.change > 0) kinds.add('cluster_distortion');
            else if (rigidReadable) kinds.add('rigid_cluster');
        }
        if (breathes) kinds.add('cluster_breathing');
        else if (significantChange(breathing) && rigidReadable) kinds.add('rigid_cluster');
        // An axis family the clusters favour is a shape distortion too.
        Object.values({ ...(local.long_axis_within_15deg || {}), ...(local.short_axis_within_15deg || {}) }).forEach((entry) => {
            if (entry.z >= CLUSTER_Z && entry.null > 0 && entry.obs / entry.null >= AXIS_ENHANCEMENT) kinds.add('cluster_distortion');
        });
        kinds.forEach((kind) => features.push({ kind, element }));
    }
    return { arm: 'A+', readings, verdict: normalizeVerdict({ features, ...symmetry.verdict }) };
};
