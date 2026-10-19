// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// The AI Copilot's analysis tools on the COMMITTED demo run (public/demo,
// GaTa4Se8 250 K), through the real engines: the PCA worker's message handler
// (sites / orientation / triplets) and the symmetry finder. What is pinned is
// what the model is told — the summary fields and the deterministic checks.

import { readFileSync } from 'node:fs';
import { describe, expect, it } from 'vitest';
import { structureFromRmc6f } from '../browserData';
import { handlePcaMessage, inProcessEngine } from '../workers/pcaKdeWorker';
import { displacementModulations } from '../workers/modulations';
import { siteDisplacementsFromRmc6f } from '../workers/pcaKde';
import modulationsFixture from './fixtures/modulations_fixture.json';
import clustersFixture from './fixtures/clusters_fixture.json';
import { CLUSTER_BECAUSE, clusterDistortion } from '../workers/clusters';
import { MODULATION_BECAUSE, SIGMA_MAX } from '../workers/modulations';
import { coordinationSummaryOf } from '../workers/triplets';
import { checkText, toolDescription } from '../llm/copilot/method';
import { createVerification } from '../llm/copilot/verification';
import { ALL_CASES, buildCase, buildConfiguration } from '../../validation/cases';
import {
    CLUSTER_SENTENCES,
    MODULATION_SENTENCES,
    buildCopilotTools,
    cartesianToUvw,
    summarizeBondAngles,
    tenDegreeBins,
    summarizeClusters,
    summarizeModulations,
    summarizeSymmetryLowering,
    summarizeOrientation,
    summarizeSites,
    summarizeSymmetry
} from '../copilotTools';

const RMC6F = readFileSync(new URL('../../public/demo/GTS_250K.rmc6f', import.meta.url), 'utf8');
const structure = structureFromRmc6f({ path: 'Demo/GTS_250K.rmc6f', text: RMC6F });
const request = inProcessEngine(RMC6F);
const elements = Object.keys(structure.elementCounts).sort();
const tools = Object.fromEntries(
    buildCopilotTools({ request, structure, elements }).map((tool) => [tool.name, tool])
);

describe('Copilot tools on the demo run', () => {
    it('offers the engine tools with the run\'s elements', () => {
        expect(Object.keys(tools)).toEqual([
            'list_sites', 'analyze_displacement_directions', 'compute_bond_angles', 'analyze_clusters',
            'analyze_modulations', 'detect_symmetry', 'measure_symmetry_lowering'
        ]);
        expect(tools.compute_bond_angles.parameters.properties.apex.enum).toEqual(['Ga', 'Se', 'Ta']);
        expect(tools.analyze_displacement_directions.parameters.properties.element.enum).toEqual(['Ga', 'Se', 'Ta', 'all']);
    });

    it('offers no engine tools without a reachable engine, and no symmetry tool without a basis', () => {
        expect(buildCopilotTools({ request: null, structure: { latticeVectors: structure.latticeVectors } })).toEqual([]);
        expect(buildCopilotTools({ request: null, structure }).map((tool) => tool.name)).toEqual(['detect_symmetry', 'measure_symmetry_lowering']);
    });

    it('list_sites: every site, ranked by |kurtosis|, and an ok check', async () => {
        const { data, checks } = await tools.list_sites.run({});
        const sites = await request('sites', {});
        expect(data.n_sites).toBe(sites.sites.length);
        expect(data.total_atoms).toBe(structure.totalAtoms);
        const first = sites.sites[0];
        expect(data.sites.find((row) => row.ref === first.referenceNumber)).toEqual(expect.objectContaining({
            site: first.element,
            count: first.count
        }));
        // Most remarkable first, so trimming to the budget drops the dullest
        // sites: split sites first (the demo has none), then by |kurtosis|.
        expect(data.sites.filter((row) => row.delta_A != null)).toEqual([]);
        const magnitudes = data.sites.map((row) => Math.abs(row.kurtosis));
        expect(magnitudes).toEqual([...magnitudes].sort((a, b) => b - a));
        expect(data.largest_abs_kurtosis_refs[0]).toBe(data.sites[0].ref);
        expect(checks).toEqual([{ level: 'ok', text: `${data.n_sites} sites, ${structure.totalAtoms} atoms.` }]);
    });

    it('list_sites: filters by element, and fails clearly for an absent one', async () => {
        const { data } = await tools.list_sites.run({ element: 'Ga' });
        expect(data.element).toBe('Ga');
        expect(data.sites.length).toBeGreaterThan(0);
        expect(data.sites.every((row) => row.site === 'Ga')).toBe(true);
        const none = summarizeSites({ sites: [{ referenceNumber: 1, element: 'Se', count: 4 }] }, { element: 'Nb' });
        expect(none.checks).toEqual([{ level: 'fail', text: 'No site holds Nb.' }]);
    });

    it('analyze_displacement_directions: one site, with crystal directions and calibrated checks', async () => {
        const { data, checks } = await tools.analyze_displacement_directions.run({ site_ref: 1 });
        expect(data.site).toBe(1);
        expect(data.used_points).toBeGreaterThan(0);
        expect(data.peak.direction_uvw).toHaveLength(3);
        expect(Math.max(...data.peak.direction_uvw.map(Math.abs))).toBe(1);
        expect(data.peak.angles_to_abc_deg).toHaveLength(3);
        expect(data.orientation_tensor.principal_axis_uvw).toHaveLength(3);
        // Exactly one verdict on the peak, calibrated against 3σ.
        const peakChecks = checks.filter((check) => /peak/i.test(check.text));
        expect(peakChecks).toHaveLength(1);
        expect(peakChecks[0].level).toBe(data.peak.significance_sigma >= 3 ? 'ok' : 'note');
    });

    it('analyze_displacement_directions: pools an element, and needs a site or an element', async () => {
        const { data } = await tools.analyze_displacement_directions.run({ element: 'Ta', frequency: 4 });
        expect(data).toMatchObject({ site: null, element: 'Ta', frequency: 4, cells: 162 });
        expect(data.used_points).toBe(structure.elementCounts.Ta);
        await expect(tools.analyze_displacement_directions.run({})).rejects.toThrow('give site_ref');
    });

    it('analyze_displacement_directions: a pooled map says its sites\' frames may cancel; one site does not', async () => {
        const pooled = await tools.analyze_displacement_directions.run({ element: 'Ta', frequency: 4 });
        const note = pooled.checks.find((check) => /^Pooled over/.test(check.text));
        expect(note.text).toMatch(/^Pooled over 16 sites\./);
        expect(note.text).toMatch(/analyze_clusters/);
        const single = await tools.analyze_displacement_directions.run({ site_ref: 1, frequency: 4 });
        expect(single.checks.some((check) => /^Pooled over/.test(check.text))).toBe(false);
    });

    it('compute_bond_angles: the Ta4 cluster gives 60° Ta–Ta–Ta angles', async () => {
        const { data, checks } = await tools.compute_bond_angles.run({
            end1: 'Ta', apex: 'Ta', end2: 'Ta', r12_min: 2.7, r12_max: 3.3
        });
        expect(data.triplet).toBe('Ta–Ta–Ta (Ta central)');
        expect(data.angle_count).toBeGreaterThan(0);
        // Most Ta atoms have their three cluster neighbours in the window (a
        // few lose one to thermal spread at 250 K): angles of ~60°.
        const counts = data.coordination.counts;
        expect(counts.indexOf(Math.max(...counts))).toBe(3);
        expect(data.coordination.mean).toBeGreaterThan(2.5);
        expect(data.modes[0].angle_deg).toBeGreaterThan(55);
        expect(data.modes[0].angle_deg).toBeLessThan(65);
        expect(data.histogram_10deg.reduce((sum, value) => sum + value, 0)).toBeCloseTo(1, 1);
        expect(checks.at(-1).level).toBe('ok');
        // The bond spread is measured, from the engine's length histogram.
        expect(data.bond_ab.sd_A).toBeGreaterThan(0.05);
        expect(data.bond_ab.sd_A).toBeLessThan(0.2);
    });

    it('analyze_clusters: the demo\'s Ta4 tetrahedra against independent atoms', async () => {
        const { data, checks } = await tools.analyze_clusters.run({ element: 'Ta', r_min: 2.65, r_max: 3.33 });
        expect(data.clusters).toMatchObject({ members: 4, edges: 6, count: 4000, templates: 4 });
        expect(data.local.ladder_A.obs).toHaveLength(6);
        expect(data.local.tetrahedral_modes.t2_trigonality.obs).toBeGreaterThan(1 / 3);
        // The RMC clusters are more rigid than independent atoms make them.
        expect(data.local.ladder_width_A.change).toBeLessThan(-0.1);
        expect(checks.find((check) => /LESS distorted/.test(check.text)).level).toBe('ok');
        expect(checks.some((check) => /compressed along <111>/.test(check.text))).toBe(true);
        // Fits, with its checks and ref, the result budget the tool asks for:
        // cut to text, its numbers could not be quoted by reference (THY-18).
        expect(JSON.stringify(data).length).toBeLessThan(3000);
        expect(JSON.stringify({ ref: 'r1', checks, data }).length).toBeLessThan(tools.analyze_clusters.resultBudget);
        expect(tools.analyze_clusters.label({ element: 'Ta', r_min: 2.65, r_max: 3.33 })).toBe('Clusters Ta · 2.65–3.33 Å');
        await expect(tools.analyze_clusters.run({ element: 'Ta', r_min: 3.3, r_max: 2.6 })).rejects.toThrow('r_min must be below r_max');
    });

    it('compute_bond_angles: a window that cuts the TaSe6 shell is refused; the whole shell is read against an octahedron', async () => {
        // gemma4:26b on Ollama picked 2.3–2.7 Å, read coordination 3.75 as
        // "incomplete octahedra", and called the answer achieved. The shell
        // has two Ta–Se peaks (2.52 and 2.74 Å): the window cuts the second.
        const cut = await tools.compute_bond_angles.run({ end1: 'Se', apex: 'Ta', end2: 'Se', r12_min: 2.3, r12_max: 2.7 });
        expect(cut.data.coordination.mean).toBeLessThan(4);
        expect(cut.checks).toHaveLength(1);
        expect(cut.checks[0].level).toBe('fail');
        expect(cut.checks[0].text).toMatch(/^Refused: Only 47 % of the Ta atoms have the most common count \(3 Se–Ta bonds\); the counts spread from 1 to 6\. The window cuts through a coordination shell\. Use the whole shell/);
        // The cut is shown; its angles are not given to quote.
        expect(cut.data.refused).toMatch(/cuts the coordination shell/);
        for (const key of ['modes', 'mean_deg', 'sd_deg', 'histogram_10deg', 'histogram_10deg_sin_corrected', 'polyhedron']) {
            expect(cut.data).not.toHaveProperty(key);
        }
        // The user may still want that window: confirmed, it is kept with a warning.
        const chosen = await tools.compute_bond_angles.run({ end1: 'Se', apex: 'Ta', end2: 'Se', r12_min: 2.3, r12_max: 2.7, window_confirmed: true });
        expect(chosen.checks[0].level).toBe('warn');
        expect(chosen.checks[0].text).toMatch(/Kept: the user chose this window/);
        expect(chosen.data.modes.length).toBeGreaterThan(0);

        const whole = await tools.compute_bond_angles.run({ end1: 'Se', apex: 'Ta', end2: 'Se', r12_min: 2.33, r12_max: 2.97 });
        expect(whole.data.coordination.counts.indexOf(Math.max(...whole.data.coordination.counts))).toBe(6);
        expect(whole.checks.map((check) => check.level)).toEqual(['ok', 'note']);
        // Read against an ideal octahedron: about a fifth of the angles are trans.
        const { cis_below_135_deg: cis, trans_above_135_deg: trans } = whole.data.polyhedron;
        expect(whole.data.polyhedron.ideal).toMatch(/octahedron/);
        expect(trans.share).toBeGreaterThan(0.15);
        expect(trans.share).toBeLessThan(0.25);
        expect(trans.mean_deg).toBeGreaterThan(150);
        expect(cis.mean_deg).toBeGreaterThan(85);
        expect(cis.mean_deg).toBeLessThan(100);
        expect(whole.checks[1].text).toMatch(/^Raw counts deplete angles near 180° \(as sin θ\)/);
        // Corrected for geometry the trans region stands out: the raw counts hide it.
        const corrected = whole.data.histogram_10deg_sin_corrected;
        expect(corrected).toHaveLength(18);
        expect(Math.max(...corrected.slice(15))).toBeGreaterThan(Math.max(...whole.data.histogram_10deg.slice(15)) * 3);
    });

    it('compute_bond_angles: an empty window fails with advice, and bad windows are refused', async () => {
        const { data, checks } = await tools.compute_bond_angles.run({
            end1: 'Ta', apex: 'Ta', end2: 'Ta', r12_min: 0.1, r12_max: 0.5
        });
        expect(data.angle_count).toBe(0);
        expect(checks).toHaveLength(1);
        expect(checks[0].level).toBe('fail');
        expect(checks[0].text).toMatch(/Widen the window/);
        await expect(tools.compute_bond_angles.run({ end1: 'Ta', apex: 'Ta', end2: 'Ta', r12_min: 3, r12_max: 2 }))
            .rejects.toThrow('r12_min must be below r12_max');
        await expect(tools.compute_bond_angles.run({ end1: 'Ta', apex: 'Ta', end2: 'Ta', r12_min: 2, r12_max: 3, r23_min: 2 }))
            .rejects.toThrow('together');
    });

    it('analyze_modulations: no cell doubling in the demo, and its box-scale waves noted', async () => {
        const { data, checks } = await tools.analyze_modulations.run({});
        expect(data.supercell).toEqual([10, 10, 10]);
        expect(data.elements.map((entry) => entry.element)).toEqual(['Ga', 'Se', 'Ta']);
        expect(data.elements[0].zone_boundary).toHaveLength(7);
        expect(checks[0].level).toBe('ok');
        expect(checks[0].text).toMatch(/^No order that doubles the cell/);
        expect(checks.some((check) => /long-wavelength wave/.test(check.text))).toBe(true);
        expect(JSON.stringify(data).length).toBeLessThan(3000);
        const one = await tools.analyze_modulations.run({ element: 'Ta' });
        expect(one.data.elements.map((entry) => entry.element)).toEqual(['Ta']);
    });

    it('tells the Bond Geometry and Displacement Directions pages what to load', () => {
        expect(tools.compute_bond_angles.pageParams({ end1: 'Se', apex: 'Ta', end2: 'Se', r12_min: 2.4, r12_max: 2.9 })).toEqual({
            end1: 'Se', apex: 'Ta', end2: 'Se', r12Min: 2.4, r12Max: 2.9, r23Min: null, r23Max: null, binWidth: 1
        });
        expect(tools.analyze_displacement_directions.pageParams({ site_ref: 23 })).toEqual({ referenceNumber: 23, frequency: 'auto' });
        expect(tools.analyze_displacement_directions.pageParams({ site_ref: 2, frequency: 8 })).toEqual({ referenceNumber: 2, frequency: 8 });
        // The page cannot pool an element: it is only opened.
        expect(tools.analyze_displacement_directions.pageParams({ element: 'Ta' })).toBeNull();
        expect(tools.list_sites.pageParams).toBeUndefined();
    });

    it('detect_symmetry: F-43m at the card default, as the Detected SG card reports', async () => {
        const { data, checks } = await tools.detect_symmetry.run({ tolerance_A: 0.2 });
        expect(data).toMatchObject({ tolerance_A: 0.2, space_group: 'F-43m', number: 216 });
        expect(data.orbits.map((orbit) => orbit.wyckoff)).toContain('4c');
        expect(checks).toEqual([
            { level: 'ok', text: 'F-43m (No. 216) holds at 0.2 Å.' },
            {
                level: 'note',
                text: 'Point group -43m: cubic; no inversion centre; achiral (not a Sohncke group); non-polar. '
                    + 'Take inversion, chirality and polarity from this, not from memory.',
                // The advice arm C's method 'off' leaves out (llm copilot/method.js).
                method: [' Take inversion, chirality and polarity from this, not from memory.']
            }
        ]);
        expect(data.point_group_facts).toBe('-43m: cubic; no inversion centre; achiral (not a Sohncke group); non-polar');
    });

    it('detect_symmetry: below the noise floor, P-42₁m — achiral and without inversion, like F-43m', async () => {
        const result = await tools.detect_symmetry.run({ tolerance_A: 0.02 });
        expect(result.data).toMatchObject({ space_group: 'P-42_1m', point_group: '-42m' });
        expect(result.facts).toEqual({
            groups: [{ symbol: 'P-42_1m', point_group: '-42m', centrosymmetric: false, sohncke: false, polar: false, crystal_system: 'tetragonal' }]
        });
    });
});

describe('Copilot summarizers', () => {
    it('cartesianToUvw keeps the sign and scales the largest component to 1', () => {
        const hexagonal = [[3, 0, 0], [-1.5, 1.5 * Math.sqrt(3), 0], [0, 0, 5]];
        expect(cartesianToUvw([0, 0, -1], hexagonal)).toEqual([0, 0, -1]);
        // Along a + b in a 120° cell: Cartesian (0.75, 1.299, 0) direction.
        const ab = [1.5, 1.5 * Math.sqrt(3), 0];
        expect(cartesianToUvw(ab.map((value) => value / 3), hexagonal)).toEqual([1, 1, 0]);
        expect(cartesianToUvw([1, 0, 0], [[1, 0, 0], [2, 0, 0], [0, 0, 1]])).toBeNull();
    });

    it('orientation checks: low statistics, no significant peak, undefined map test, asymmetry', () => {
        const { checks } = summarizeOrientation({
            peakDirection: [1, 0, 0], peakSignificance: 1.2, peakEnhancement: 1.4, usedPoints: 40,
            mapSignificance: null, antipodalAsymmetrySignificant: true, antipodalAsymmetryZ: 3.4,
            parseWarning: '2 lines skipped'
        }, null);
        expect(checks.map((check) => check.level)).toEqual(['warn', 'warn', 'note', 'note', 'note']);
        expect(checks[2].text).toBe('No significant peak direction (1.2σ < 3σ).');
    });

    it('bond-angle checks: a window that misses most central atoms is refused', () => {
        const { checks } = summarizeBondAngles({
            triplet: ['Se', 'Ta', 'Se'], angleCount: 2, apexCount: 10, coordination: [8, 0, 2], coordinationSummary: coordinationSummaryOf([8, 0, 2]),
            binCenters: [89.5, 90.5], counts: [1, 1], meanAngle: 90, stdAngle: 0.5,
            lengths12: { meanLength: 2.6, uniqueBonds: 4 }, sharedEnds: true
        });
        expect(checks).toEqual([{
            level: 'fail',
            text: 'Refused: 80 % of the Ta atoms have no Se–Ta bond in the window: it misses most of the shell. '
                + 'Use the whole shell (find_gr_peaks\' first_shell.window_A). If the user wants this window, ask them (ask_user) '
                + 'and call again with window_confirmed: true.',
            // The advice arm C's method 'off' leaves out (llm copilot/method.js).
            method: [' Use the whole shell (find_gr_peaks\' first_shell.window_A).']
        }]);
    });

    it('bond-angle checks read the engine\'s coordinationSummary through the shell module\'s counts rule, and fail loudly without it', () => {
        const payload = {
            triplet: ['Se', 'Ta', 'Se'], angleCount: 20, apexCount: 10, coordination: [0, 0, 4, 3, 3],
            binCenters: [89.5, 90.5], counts: [10, 10], meanAngle: 90, stdAngle: 0.5,
            lengths12: { meanLength: 2.6, uniqueBonds: 14 }, sharedEnds: true
        };
        // 4 of 10 share the modal count (2): below 60 %, refused; the tie-free summary is the engine's.
        const refusedResult = summarizeBondAngles({ ...payload, coordinationSummary: coordinationSummaryOf(payload.coordination) });
        expect(refusedResult.data.refused).toMatch(/cuts the coordination shell/);
        expect(refusedResult.checks[0].text).toMatch(/^Refused: Only 40 % of the Ta atoms have the most common count \(2 Se–Ta bonds\); the counts spread from 2 to 4\./);
        expect(refusedResult.data.coordination.mean).toBe(2.9);
        expect(() => summarizeBondAngles(payload)).toThrow(/coordinationSummary/);
    });

    it('bond-angle reference: only an A–B–A triplet in one window is read against a polyhedron (review 2026-10-12 #5)', () => {
        // Six Se per Ta, every angle near 90°: a whole octahedral shell's
        // numbers, but for a mixed or a split triplet the angles are a subset.
        const base = {
            angleCount: 15, apexCount: 1, coordination: [0, 0, 0, 0, 0, 0, 1], coordinationSummary: coordinationSummaryOf([0, 0, 0, 0, 0, 0, 1]),
            binCenters: [89.5, 90.5], counts: [8, 7], meanAngle: 90, stdAngle: 0.5,
            lengths12: { meanLength: 2.6, uniqueBonds: 6 }
        };
        const mixed = summarizeBondAngles({ ...base, triplet: ['Se', 'Ta', 'Te'], sharedEnds: false, bond23: [2.4, 2.9] });
        expect(mixed.data).not.toHaveProperty('polyhedron');
        expect(mixed.checks.map((check) => check.level)).toEqual(['ok']);
        const split = summarizeBondAngles({ ...base, triplet: ['Se', 'Ta', 'Se'], sharedEnds: false, bond12: [2.3, 2.6], bond23: [2.6, 2.9] });
        expect(split.data).not.toHaveProperty('polyhedron');
        // One shell: read against the octahedron; no trans angles is "none", never "null°".
        const shell = summarizeBondAngles({ ...base, triplet: ['Se', 'Ta', 'Se'], sharedEnds: true });
        expect(shell.data.polyhedron.trans_above_135_deg).toEqual({ share: 0, mean_deg: null });
        expect(shell.checks[1].text).toMatch(/0 % of these angles lie above 135° \(mean none\)/);
        expect(shell.checks[1].text).not.toMatch(/null/);
    });

    it('bond-angle histograms: wide engine bins are shared by overlap, not combed (review 2026-10-12 #5)', () => {
        // A flat distribution in 15° bins: each 10° bin gets two thirds of one.
        const centers = Array.from({ length: 12 }, (_, i) => 7.5 + 15 * i);
        const flat = tenDegreeBins(centers, centers.map(() => 150));
        flat.forEach((value) => expect(value).toBeCloseTo(100, 9));
        // 30° bins too, and 1° bins as before (whole bins).
        const wide = tenDegreeBins([15, 45, 75, 105, 135, 165], [300, 300, 300, 300, 300, 300]);
        wide.forEach((value) => expect(value).toBeCloseTo(100, 9));
        const fine = Array.from({ length: 180 }, (_, i) => i + 0.5);
        expect(tenDegreeBins(fine, fine.map((_, i) => i))).toEqual(
            Array.from({ length: 18 }, (_, j) => Array.from({ length: 10 }, (_, k) => 10 * j + k).reduce((sum, value) => sum + value, 0))
        );
        const summary = summarizeBondAngles({
            triplet: ['Se', 'Ta', 'Se'], sharedEnds: true, angleCount: 1800, apexCount: 100, coordination: [0, 0, 0, 0, 0, 0, 100],
            coordinationSummary: coordinationSummaryOf([0, 0, 0, 0, 0, 0, 100]),
            binCenters: centers, counts: centers.map(() => 150), meanAngle: 90, stdAngle: 50, lengths12: { meanLength: 2.6, uniqueBonds: 600 }
        });
        expect(summary.data.histogram_10deg.every((value) => value > 0.05)).toBe(true);
    });

    it('cluster checks: an extended network and an empty window fail with advice, and give no facts', () => {
        const parsed = siteDisplacementsFromRmc6f(clustersFixture.cases.find((entry) => entry.name === 'trigonal').text);
        const extended = summarizeClusters(clusterDistortion(parsed, 'Se', [7.5, 8]));
        expect(extended.checks).toHaveLength(1);
        expect(extended.checks[0].level).toBe('fail');
        expect(extended.checks[0].text).toMatch(/extended network/);
        expect(extended.facts).toBeUndefined();
        const empty = summarizeClusters(clusterDistortion(parsed, 'Nb', [0.5, 1]));
        expect(empty.data.sites_without_bond).toBe(4);
        expect(empty.checks[0].text).toMatch(/^No two Nb sites/);
        expect(empty.facts).toBeUndefined();
    });

    it('the sentence tables: every because the engines can give has an entry (a check, or null when another finding says it)', () => {
        Object.entries(CLUSTER_BECAUSE).forEach(([name, codes]) => {
            expect(Object.keys(CLUSTER_SENTENCES[name] || {}).sort(), name).toEqual([...codes].sort());
        });
        expect(Object.keys(CLUSTER_SENTENCES).sort()).toEqual(Object.keys(CLUSTER_BECAUSE).sort());
        // A point is said as a row, not a sentence.
        Object.entries(MODULATION_BECAUSE).filter(([name]) => name !== 'point').forEach(([name, codes]) => {
            expect(Object.keys(MODULATION_SENTENCES[name] || {}).sort(), name).toEqual([...codes].sort());
        });
        expect(Object.keys(MODULATION_SENTENCES).sort()).toEqual(['cellDoubling', 'generalOrder', 'longWave']);
        // A because the table does not know is a broken contract with the engine: it throws, never goes unsaid.
        const payload = clusterDistortion(siteDisplacementsFromRmc6f(clustersFixture.cases.find((entry) => entry.name === 'trigonal').text), 'Nb', [2.5, 3.5]);
        const broken = { ...payload, verdict: { ...payload.verdict, breathing: { ...payload.verdict.breathing, because: 'invented_code' } } };
        expect(() => summarizeClusters(broken)).toThrow(/breathing finding gave 'invented_code', which has no sentence/);
    });

    it('cluster facts are the engine\'s answers', () => {
        const payload = clusterDistortion(siteDisplacementsFromRmc6f(RMC6F), 'Ta', [2.65, 3.33]);
        const { facts } = summarizeClusters(payload);
        const { verdict } = payload;
        expect(facts.cluster).toEqual({
            element: 'Ta',
            members: 4,
            window: [2.65, 3.33],
            distorted: verdict.distorted.answer,
            static: verdict.staticDistortion.answer,
            breathing: verdict.breathing.answer,
            rigid: verdict.rigid.answer,
            split: verdict.splitMemberSites.answer,
            families: ['elongated <111>', 'compressed <111>']
        });
        expect(facts.cluster).toMatchObject({ distorted: true, static: true, breathing: true, rigid: false, split: false });
    });

    it('modulation sentences: a zone-boundary point the engine calls ordered (σ = +∞ reads > 8σ)', () => {
        const bar = { ratio: 3, sigma: 3 };
        const point = (q, ratio, percentile, sigma, ordered, because, kind = 'zone_boundary') => ({ q, ratio, percentile, rank: 1, kind, sigma, ordered, because, bar });
        const ordered = point([0, 0, 0.5], 204, 99.53, SIGMA_MAX, true, 'ordered');
        const finding = (answer, because, decidedBy = [], nearest = null) => ({ answer, because, decidedBy, nearest });
        const { data, checks, facts } = summarizeModulations({
            supercell: [6, 6, 6],
            wavevectors: 215,
            untestable: [],
            family: { zoneBoundaryPoints: 2, comparedModes: 98 },
            rules: { ratio: 3, sigma: 3, longWaveRatio: 5, longWaveQ: 0.25, sigmaMax: SIGMA_MAX },
            elements: [{
                element: 'Nb',
                sites: 4,
                medianPower: 0.85,
                comparedModes: 98,
                zoneBoundary: [ordered, point([0.5, 0, 0], 0.9, 40, -2.4, false, 'below_ratio')],
                strongest: [{ ...ordered, sigma: null, ordered: null, because: 'zone_boundary', bar: null }, point([0.1, 0, 0], 2, 90, null, null, 'long_wave', 'long_wave')],
                verdict: {
                    cellDoubling: finding(true, 'ordered_point', [[0, 0, 0.5]]),
                    generalOrder: finding(null, 'no_general_mode'),
                    longWave: finding(false, 'strongest_not_long_wave')
                }
            }]
        });
        expect(checks).toHaveLength(1);
        expect(checks[0]).toEqual({
            level: 'ok',
            text: 'Ordered Nb displacements at q = (0, 0, ½): 204× the median power, > 8σ against independent displacements '
                + '(corrected over 2 zone-boundary points) — a superstructure doubling the cell along c.'
        });
        // Ratios to two significant digits (204 reads 200), σ past SIGMA_MAX as '>8'.
        expect(data.elements[0].zone_boundary).toEqual([['(0, 0, ½)', 200, 99.5, '>8'], ['(½, 0, 0)', 0.9, 40, -2.4]]);
        expect(data.elements[0].strongest.map((row) => row[3])).toEqual([null, null]);
        expect(facts.modulation).toEqual({ elements: ['Nb'], ordered: [{ element: 'Nb', q: [0, 0, 0.5] }], general: [] });
        const odd = summarizeModulations({
            supercell: [3, 4, 4], wavevectors: 47, untestable: [[1, 0, 0]], family: { zoneBoundaryPoints: 0, comparedModes: 0 },
            rules: { ratio: 3, sigma: 3, longWaveRatio: 5, longWaveQ: 0.25, sigmaMax: SIGMA_MAX },
            elements: [{ element: 'Nb', sites: 0, verdict: { cellDoubling: finding(null, 'no_displacements'), generalOrder: finding(null, 'no_displacements'), longWave: finding(null, 'no_displacements') } }]
        });
        expect(odd.checks[0].level).toBe('fail');
    });

    it('modulation checks: order under a stronger long wave is still found (review 2026-10-12 #6)', () => {
        // Ti order at (0, 0, ½), 55× the median, under a long wave at (⅙, 0, 0)
        // and its −q partner: percentile 98.6 of all q, so it was missed.
        const { text } = modulationsFixture.cases.find((entry) => entry.name.startsWith('long-wave'));
        const payload = displacementModulations(siteDisplacementsFromRmc6f(text), 'Ti');
        const { checks } = summarizeModulations(payload);
        expect(checks[0].level).toBe('ok');
        expect(checks[0].text).toMatch(/^Ordered Ti displacements at q = \(0, 0, ½\): \d+(\.\d+)?× the median power, (> 8|\d+(\.\d+)?)σ against independent displacements/);
        expect(checks.some((check) => /long-wavelength wave/.test(check.text))).toBe(true);
    });

    it('modulation sentences: the planted order said, and only it; the decoy said as a long wave (dev-zone-boundary); none where there is none', () => {
        const planted = buildCase(ALL_CASES.find((item) => item.id === 'dev-zone-boundary'), { partials: false });
        const { checks, facts } = summarizeModulations(displacementModulations(siteDisplacementsFromRmc6f(planted.rmc6f)));
        expect(facts.modulation.ordered).toEqual([{ element: 'Ti', q: [0, 0, 0.5] }]);
        expect(checks[0].text).toMatch(/^Ordered Ti displacements at q = \(0, 0, ½\)/);
        expect(checks.filter((check) => /^Ordered/.test(check.text))).toHaveLength(1);
        expect(checks.some((check) => /Sr's strongest mode is a long-wavelength wave/.test(check.text))).toBe(true);
        const thermal = buildCase(ALL_CASES.find((item) => item.id === 'dev-thermal-null'), { partials: false });
        const none = summarizeModulations(displacementModulations(siteDisplacementsFromRmc6f(thermal.rmc6f)));
        expect(none.checks.find((check) => check.level === 'ok').text).toMatch(/^No order that doubles the cell: no zone-boundary point stands out of the diffuse background \(the strongest, /);
        expect(none.facts.modulation).toMatchObject({ ordered: [], general: [] });
    });

    it('displacement directions: a thermal ellipsoid is anisotropic, not a preferred direction (review L3)', async () => {
        // An O flattened across its bond, thermal motion only: the map test
        // passes (7.2σ), the peak test does not (2.3σ). The direction fact
        // used to count the map test, so "preferred direction" went unflagged.
        const built = buildCase(ALL_CASES.find((item) => item.id === 'dev-thermal-null'), { partials: false });
        const payload = await handlePcaMessage({ kind: 'orientation', referenceNumber: 3, geometry: false }, async () => built.rmc6f);
        const { checks, facts } = summarizeOrientation(payload, null);
        expect(payload.peakSignificance).toBeLessThan(3);
        expect(payload.mapSignificance).toBeGreaterThan(3);
        expect(facts.direction).toMatchObject({ subject: 'site 3', significant: false, anisotropic: true });
        expect(checks.find((check) => /anisotropic/.test(check.text)).text)
            .toMatch(/^The displacements are anisotropic: the map departs from isotropy at 7\.\d+σ\. Thermal ellipsoids pass this test too: it is not a preferred direction/);
        // The answer check on it (the verification: the step recorded with its facts, then the answer judged).
        const verification = createVerification({ offered: [] });
        verification.record({ name: 'analyze_displacement_directions', args: { site_ref: 3 }, status: 'ok', error: null, executed: true, facts }, () => ({ data: null }));
        const [finding] = verification.judge('Site 3 has a preferred displacement direction along [1 0 0].').done.audit;
        expect(finding.rule).toBe('direction');
        expect(finding.text).toMatch(/found a peak at ≥ 3σ \(site 3: peak 2\.\d+σ, map 7\.\d+σ\)\. A significant map test is anisotropy/);
        // A significant peak still counts, with or without the map.
        expect(summarizeOrientation({ ...payload, peakSignificance: 4.1 }, null).facts.direction.significant).toBe(true);
    });

    it('cluster summary: "LESS distorted" does not deny a distortion the axis families show (Jahn–Teller octahedra)', () => {
        const built = buildCase(ALL_CASES.find((item) => item.id === 'dev-jahn-teller'), { partials: false });
        const { checks } = summarizeClusters(clusterDistortion(siteDisplacementsFromRmc6f(built.rmc6f), 'F', [2.5, 3.25]));
        const less = checks.find((check) => /LESS distorted/.test(check.text));
        expect(less.text).toMatch(/yet more clusters than by chance are (elongated|compressed) along <100>/);
        expect(less.text).not.toMatch(/no correlated shape distortion/);
        expect(checks.some((check) => /^More clusters are (elongated|compressed) along <100>/.test(check.text))).toBe(true);
    });

    it('cluster summary: breathing is size, not shape; split sites are named (review L2)', async () => {
        // Mo3 triangles whose size fluctuates collectively (sd 4 %), shape
        // regular: the engine reads the shape at each cluster's own size
        // (its tests pin the decision), and the tool says so.
        const breathing = buildCase(ALL_CASES.find((item) => item.id === 'dev-breathing'), { partials: false });
        const ask = (built) => inProcessEngine(built.rmc6f);
        const triangles = buildCopilotTools({ request: ask(breathing), structure: null, elements: ['Mo', 'O'] })
            .find((tool) => tool.name === 'analyze_clusters');
        const { data, checks, facts } = await triangles.run({ element: 'Mo', r_min: 2.4, r_max: 3.21 });
        expect(Math.abs(data.local.ladder_width_A.z)).toBeLessThan(3);
        expect(checks.find((check) => check.level === 'ok').text).toMatch(/^No cluster shape distortion beyond independent atoms \(at each cluster's own size\)/);
        expect(checks.some((check) => /^The clusters breathe: their size \(mean edge\) fluctuates \d+ % more/.test(check.text))).toBe(true);
        expect(checks.some((check) => /LESS distorted|More clusters are/.test(check.text))).toBe(false);
        expect(facts.cluster).toMatchObject({ element: 'Mo', distorted: false, breathing: true, rigid: false, split: false });

        // Jahn–Teller F6 octahedra: the F sites have two positions each (the
        // engine finds them from its own clouds, with no second request), so
        // the null holds the distortion's amplitude; the pattern shows it.
        const jt = buildCase(ALL_CASES.find((item) => item.id === 'dev-jahn-teller'), { partials: false });
        const requests = [];
        const counting = (kind, params) => {
            requests.push(kind);
            return ask(jt)(kind, params);
        };
        const octahedra = buildCopilotTools({ request: counting, structure: null, elements: ['Cu', 'F', 'K'] })
            .find((tool) => tool.name === 'analyze_clusters');
        const split = await octahedra.run({ element: 'F', r_min: 2.5, r_max: 3.25 });
        expect(requests).toEqual(['clusters']);
        expect(split.data.split_member_sites.length).toBeGreaterThanOrEqual(3);
        expect(split.checks.find((check) => /^Split sites:/.test(check.text)).text).toMatch(/already holds the distortion's amplitude/);
        expect(split.checks.some((check) => /the split sites' positions combine in a pattern that keeps the size/.test(check.text))).toBe(true);
        expect(split.facts.cluster).toMatchObject({ distorted: true, split: true, rigid: false });
    });

    it('measure_symmetry_lowering: the coherent shift the 0.2 Å card hides, and none where there is none (review L4)', async () => {
        const tool = (id) => {
            const built = buildCase(ALL_CASES.find((item) => item.id === id), { partials: false });
            const parsed = structureFromRmc6f({ path: `${id}.rmc6f`, text: built.rmc6f });
            return buildCopilotTools({ request: null, structure: parsed, elements: [] }).find((entry) => entry.name === 'measure_symmetry_lowering');
        };
        // Every Ta shifted 0.07 Å along [001] (P4mm): at 0.2 Å the card says Pm-3m.
        const lowering = await tool('dev-symmetry-lowering').run({});
        expect(lowering.data.parent).toMatchObject({ tolerance_A: 0.2, space_group: 'Pm-3m' });
        expect(lowering.data.lower.space_group).toBe('P4mm');
        expect(lowering.data.sites[0]).toMatchObject({ element: 'Ta', above_noise_floor: true, nearest_low_index: { uvw: [0, 0, 1] } });
        expect(Math.abs(lowering.data.sites[0].shift_A - 0.07)).toBeLessThan(0.01);
        expect(lowering.data.sites.slice(1).every((site) => !site.above_noise_floor)).toBe(true);
        expect(lowering.checks[0].text).toMatch(/^Ta shifted 0\.0\d+ Å along \[0 0 1\] \([\d.]+° off\) from the Pm-3m positions .* \(the lower group P4mm at 0\.0\d+ Å\)/);
        expect(lowering.checks[1].text).toBe('Pm-3m → P4mm loses the inversion centre and becomes polar.');
        expect(lowering.facts.lowering).toMatchObject({ parent: 'Pm-3m', lower: 'P4mm', shifted: [{ element: 'Ta', direction_uvw: [0, 0, 1] }] });
        // Thermal motion only: nothing above the noise floor.
        const none = await tool('dev-thermal-null').run({});
        expect(none.facts.lowering.shifted).toEqual([]);
        expect(none.checks[0].text).toMatch(/^No site is shifted from the Pm-3m positions above the 0\.0\d+ Å noise floor/);
        expect(summarizeSymmetryLowering(null, null).checks[0].level).toBe('fail');
    });

    it('modulation sentences: correlation between cells at a general wavevector (review L4)', () => {
        // q = (3/8, 0, 0) on an 8-cell axis: neither a zone-boundary point nor a long wave.
        const config = buildConfiguration({
            cell: [4, 4, 4], supercell: [8, 8, 8], seed: 11,
            sites: [{ el: 'Ti', frac: [0.5, 0.5, 0.5] }, { el: 'O', frac: [0.5, 0.5, 0] }],
            shift: (index, _, rng, n) => {
                const thermal = [0, 1, 2].map(() => 0.06 * Math.sqrt(-2 * Math.log(rng() || 1e-12)) * Math.cos(2 * Math.PI * rng()));
                return index === 0 ? [thermal[0], thermal[1], thermal[2] + 0.08 * Math.cos(2 * Math.PI * 0.375 * n[0])] : thermal;
            }
        });
        const { checks, facts, data } = summarizeModulations(displacementModulations(siteDisplacementsFromRmc6f(config.text)));
        expect(facts.modulation.ordered).toEqual([]);
        expect(facts.modulation.general).toEqual([{ element: 'Ti', q: [expect.closeTo(0.375, 6), 0, 0] }]);
        expect(checks.find((check) => /^Correlated Ti displacements between cells at q = \(0\.375, 0, 0\)/.test(check.text)).text)
            .toMatch(/against independent displacements \(corrected over \d+ compared modes\) — order at a general wavevector .* not a cell doubling\.$/);
        expect(data.elements.find((entry) => entry.element === 'Ti').strongest[0][3]).not.toBeNull();
    });

    it('cluster summary: other_sizes never lists the analysed kind, on a tie too (review 2026-10-12 #17)', () => {
        const { text } = clustersFixture.cases.find((entry) => entry.name === 'tie');
        const payload = clusterDistortion(siteDisplacementsFromRmc6f(text), 'Mo', [1.8, 2.2]);
        expect(payload.kinds[0]).toMatchObject({ members: 2, clusters: 64 });
        const { data } = summarizeClusters(payload);
        expect(data.clusters).toMatchObject({ members: 2, count: 64 });
        expect(data.clusters.other_sizes).toEqual([{ members: 3, count: 64 }]);
    });

    it('symmetry checks: not analysed, undetermined, lower bound', () => {
        expect(summarizeSymmetry({ skipped: true, reason: 'too many', spaceGroup: 'not analysed', orbits: [] }, 0.2).checks[0])
            .toEqual({ level: 'warn', text: 'Not analysed: too many' });
        expect(summarizeSymmetry({ spaceGroup: 'undetermined', orbits: [] }, 0.2).checks[0].level).toBe('note');
        expect(summarizeSymmetry({ spaceGroup: '≥ -43m', orbits: [] }, 0.2).checks[0].text).toMatch(/lower bound/);
        expect(summarizeSymmetry(null, 0.2).checks[0].level).toBe('fail');
    });
});

// The method advice in the engine tools (llm copilot/method.js, plan item
// 1.8): the validation's arm C (method 'off') reads each description without
// it and each check without it — a refusal still refuses, a fact stays. With
// method 'full' (the dock, arm B) every text is as before; the model-visible
// snapshots pin that, and arm C's tools (harness/tools.arm-C.json).
describe('the engine tools\' method advice: what arm C reads without it', () => {
    // Every tool the Copilot offers a model, by name.
    const TOOL_NAME = /\b(?:read_run_data|assess_convergence|find_gr_peaks|compare_fit|list_sites|analyze_[a-z_]+|compute_bond_angles|detect_symmetry|measure_symmetry_lowering|ask_user|read_skill)\b/g;
    const toolNames = (text) => text.match(TOOL_NAME) || [];

    it('each description keeps what the tool does and returns: no other tool named but site_ref\'s source', () => {
        const advised = Object.values(tools).filter((tool) => tool.method != null).map((tool) => tool.name);
        expect(advised).toEqual(['list_sites', 'compute_bond_angles', 'analyze_clusters', 'analyze_modulations', 'detect_symmetry', 'measure_symmetry_lowering']);
        for (const tool of Object.values(tools)) {
            const off = toolDescription(tool, 'off');
            expect(off).toBe(tool.description);
            // site_ref's values are list_sites' refs: the parameter's meaning, not advice.
            expect(toolNames(off)).toEqual(tool.name === 'analyze_displacement_directions' ? ['list_sites'] : []);
        }
        expect(toolDescription(tools.compute_bond_angles, 'full')).toMatch(/tetrahedron\. Use the whole first coordination shell \(find_gr_peaks' first_shell\.window_A\): /);
        expect(toolDescription(tools.analyze_clusters, 'full'))
            .toContain('in the average structure (take the window from find_gr_peaks: first_shell.window_A of that element pair). Every cluster');
        expect(toolDescription(tools.detect_symmetry, 'full')).toMatch(/Wyckoff orbits\. Compare tolerances to see which symmetry the local distortion breaks\. Does not change the card\.$/);
        // The refusal stays stated in arm C: on window_confirmed.
        expect(tools.compute_bond_angles.parameters.properties.window_confirmed.description).toMatch(/A window that cuts the coordination shell is refused unless it is set\./);
    });

    it('compute_bond_angles: the cut-window refusal still refuses; only which window to use goes', async () => {
        const cut = await tools.compute_bond_angles.run({ end1: 'Se', apex: 'Ta', end2: 'Se', r12_min: 2.3, r12_max: 2.7 });
        const [check] = cut.checks;
        expect(checkText(check, 'full')).toBe(check.text);
        expect(checkText(check, 'off')).toBe('Refused: Only 47 % of the Ta atoms have the most common count (3 Se–Ta bonds); the counts spread from 1 to 6. '
            + 'The window cuts through a coordination shell. If the user wants this window, ask them (ask_user) and call again with window_confirmed: true.');
        expect(cut.data.refused).toMatch(/cuts the coordination shell/);
    });

    it('compute_bond_angles: no angles — the failure and the plain remedy, not where the window comes from', async () => {
        const { checks } = await tools.compute_bond_angles.run({ end1: 'Ta', apex: 'Ta', end2: 'Ta', r12_min: 0.1, r12_max: 0.5 });
        expect(checks[0].text).toBe('No angles: no central atom has two neighbours in the window(s). Widen the window '
            + '(find_gr_peaks gives the whole first shell, first_shell.window_A) or check the elements.');
        expect(checkText(checks[0], 'off')).toBe('No angles: no central atom has two neighbours in the window(s). Widen the window or check the elements.');
    });

    it('analyze_displacement_directions: a pooled map keeps its caveat; which analysis to run instead goes', async () => {
        const pooled = await tools.analyze_displacement_directions.run({ element: 'Ta', frequency: 4 });
        const note = pooled.checks.find((check) => /^Pooled over/.test(check.text));
        expect(note.text).toMatch(/cluster distortion\. Test single sites \(site_ref\), or analyze_clusters for a cluster's shape\.$/);
        expect(checkText(note, 'off')).toBe('Pooled over 16 sites. Where their local axes differ (the four vertices of a tetrahedral cluster), a direction '
            + 'each site prefers in its own frame cancels here, so an isotropic pooled map is not evidence against a per-site or cluster distortion.');
    });

    it('the reading hints that name no tool are method too: arm C reads the point group\'s facts and the octahedron counts alone', async () => {
        const symmetry = await tools.detect_symmetry.run({ tolerance_A: 0.2 });
        const pointGroup = symmetry.checks.find((check) => /^Point group/.test(check.text));
        expect(checkText(pointGroup, 'full')).toBe('Point group -43m: cubic; no inversion centre; achiral (not a Sohncke group); non-polar. '
            + 'Take inversion, chirality and polarity from this, not from memory.');
        expect(checkText(pointGroup, 'off')).toBe('Point group -43m: cubic; no inversion centre; achiral (not a Sohncke group); non-polar.');
        const whole = await tools.compute_bond_angles.run({ end1: 'Se', apex: 'Ta', end2: 'Se', r12_min: 2.33, r12_max: 2.97 });
        const octahedron = whole.checks.find((check) => /^Raw counts deplete/.test(check.text));
        expect(octahedron.text).toMatch(/at 90°\. Read the shape from histogram_10deg_sin_corrected \(1 = random directions\), not from the raw modes\.$/);
        expect(checkText(octahedron, 'off')).toMatch(/^Raw counts deplete angles near 180° \(as sin θ\): .* against an ideal octahedron's 20 % trans at 180° and 80 % cis at 90°\.$/);
    });

    it('analyze_clusters: the window\'s failures and the split sites\' caveat, without which window or analysis', () => {
        const extended = CLUSTER_SENTENCES['window.extendedNetwork'].extended_network({ element: 'O', lo: 2.4, hi: 3.13 });
        expect(extended.text).toBe('The 2.4–3.13 Å window links the O atoms into an extended network (chains or a 3-D net), not finite clusters: '
            + 'narrow it to the intra-cluster bond (find_gr_peaks: first_shell.window_A of O-O).');
        expect(checkText(extended, 'off')).toBe('The 2.4–3.13 Å window links the O atoms into an extended network (chains or a 3-D net), not finite clusters: narrow it.');
        const none = CLUSTER_SENTENCES['window.extendedNetwork'].no_clusters({ element: 'Ta', lo: 0.1, hi: 0.5 });
        expect(none.text).toBe('No two Ta sites of the average structure lie 0.1–0.5 Å apart: no clusters. Widen the window (find_gr_peaks gives the Ta-Ta first shell).');
        expect(checkText(none, 'off')).toBe('No two Ta sites of the average structure lie 0.1–0.5 Å apart: no clusters. Widen the window.');
        const split = CLUSTER_SENTENCES.splitMemberSites.half_or_more_split({ verdict: { splitMemberSites: { sites: [3, 5], memberSites: [3, 5, 7] } } });
        expect(split.text).toMatch(/does not make the clusters regular\. Read the pattern \(axis families, trigonality\), the partial g\(r\) shell \(find_gr_peaks\) and the bond angles\.$/);
        expect(checkText(split, 'off')).toMatch(/^Split sites: 2 of the 3 member sites \(3, 5\) have two positions .* does not make the clusters regular\.$/);
        for (const check of [extended, none, split]) expect(toolNames(checkText(check, 'off'))).toEqual([]);
    });
});
