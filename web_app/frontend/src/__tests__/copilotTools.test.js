// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// The AI Copilot's analysis tools on the COMMITTED demo run (public/demo,
// GaTa4Se8 250 K), through the real engines: the PCA worker's message handler
// (sites / orientation / triplets) and the symmetry finder. What is pinned is
// what the model is told — the summary fields and the deterministic checks.

import { readFileSync } from 'node:fs';
import { describe, expect, it } from 'vitest';
import { structureFromRmc6f } from '../browserData';
import { handlePcaMessage } from '../workers/pcaKdeWorker';
import {
    buildCopilotTools,
    cartesianToUvw,
    summarizeBondAngles,
    summarizeClusters,
    summarizeModulations,
    summarizeOrientation,
    summarizeSites,
    summarizeSymmetry
} from '../copilotTools';

const RMC6F = readFileSync(new URL('../../public/demo/GTS_250K.rmc6f', import.meta.url), 'utf8');
const structure = structureFromRmc6f({ path: 'Demo/GTS_250K.rmc6f', text: RMC6F });
const request = (kind, params) => handlePcaMessage({ kind, ...params }, async () => RMC6F);
const elements = Object.keys(structure.elementCounts).sort();
const tools = Object.fromEntries(
    buildCopilotTools({ request, structure, elements }).map((tool) => [tool.name, tool])
);

describe('Copilot tools on the demo run', () => {
    it('offers the engine tools with the run\'s elements', () => {
        expect(Object.keys(tools)).toEqual([
            'list_sites', 'analyze_displacement_directions', 'compute_bond_angles', 'analyze_clusters',
            'analyze_modulations', 'detect_symmetry'
        ]);
        expect(tools.compute_bond_angles.parameters.properties.apex.enum).toEqual(['Ga', 'Se', 'Ta']);
        expect(tools.analyze_displacement_directions.parameters.properties.element.enum).toEqual(['Ga', 'Se', 'Ta', 'all']);
    });

    it('offers no engine tools without a reachable engine, and no symmetry tool without a basis', () => {
        expect(buildCopilotTools({ request: null, structure: { latticeVectors: structure.latticeVectors } })).toEqual([]);
        expect(buildCopilotTools({ request: null, structure }).map((tool) => tool.name)).toEqual(['detect_symmetry']);
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
        // Most remarkable first, so trimming to the budget drops the dullest sites.
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
        // Fits the per-result budget the agent loop gives a tool.
        expect(JSON.stringify(data).length).toBeLessThan(3000);
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
                    + 'Take inversion, chirality and polarity from this, not from memory.'
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
            triplet: ['Se', 'Ta', 'Se'], angleCount: 2, apexCount: 10, coordination: [8, 0, 2],
            binCenters: [89.5, 90.5], counts: [1, 1], meanAngle: 90, stdAngle: 0.5,
            lengths12: { meanLength: 2.6, uniqueBonds: 4 }, sharedEnds: true
        });
        expect(checks).toEqual([{
            level: 'fail',
            text: 'Refused: 80 % of the Ta atoms have no Se–Ta bond in the window: it misses most of the shell. '
                + 'Use the whole shell (find_gr_peaks\' first_shell.window_A). If the user wants this window, ask them (ask_user) '
                + 'and call again with window_confirmed: true.'
        }]);
    });

    it('cluster checks: an extended network and an empty window fail with advice', () => {
        const extended = summarizeClusters({ element: 'Se', window: [7.5, 8], sites: [7], extended: true, templates: [], kinds: [] });
        expect(extended.checks).toHaveLength(1);
        expect(extended.checks[0].level).toBe('fail');
        expect(extended.checks[0].text).toMatch(/extended network/);
        const empty = summarizeClusters({ element: 'Nb', window: [0.5, 1], sites: [1, 2], isolatedSites: 2, extended: false, templates: [], kinds: [] });
        expect(empty.data.sites_without_bond).toBe(2);
        expect(empty.checks[0].text).toMatch(/^No two Nb sites/);
    });

    it('modulation checks: a zone-boundary point standing out of the background is order', () => {
        const point = (q, ratio, percentile) => ({ q, ratio, percentile });
        const { checks } = summarizeModulations({
            supercell: [6, 6, 6],
            wavevectors: 215,
            untestable: [],
            elements: [{
                element: 'Nb',
                sites: 4,
                zoneBoundary: [point([0, 0, 0.5], 204, 99.53), point([0.5, 0, 0], 0.9, 40)],
                strongest: [point([0, 0, 0.5], 204, 99.53), point([0.1, 0, 0], 2, 90)]
            }]
        });
        expect(checks).toHaveLength(1);
        expect(checks[0]).toEqual({
            level: 'ok',
            text: 'Ordered Nb displacements at q = (0, 0, ½): 204× the median power, above 99.5 % of all wavevectors — a superstructure doubling the cell along c.'
        });
        const odd = summarizeModulations({ supercell: [3, 4, 4], wavevectors: 47, untestable: [[1, 0, 0]], elements: [{ element: 'Nb', sites: 0 }] });
        expect(odd.checks[0].level).toBe('fail');
    });

    it('symmetry checks: not analysed, undetermined, lower bound', () => {
        expect(summarizeSymmetry({ skipped: true, reason: 'too many', spaceGroup: 'not analysed', orbits: [] }, 0.2).checks[0])
            .toEqual({ level: 'warn', text: 'Not analysed: too many' });
        expect(summarizeSymmetry({ spaceGroup: 'undetermined', orbits: [] }, 0.2).checks[0].level).toBe('note');
        expect(summarizeSymmetry({ spaceGroup: '≥ -43m', orbits: [] }, 0.2).checks[0].text).toMatch(/lower bound/);
        expect(summarizeSymmetry(null, 0.2).checks[0].level).toBe('fail');
    });
});
