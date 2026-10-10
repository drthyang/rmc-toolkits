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
        expect(Object.keys(tools)).toEqual(['list_sites', 'analyze_displacement_directions', 'compute_bond_angles', 'detect_symmetry']);
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
        expect(checks).toEqual([{ level: 'ok', text: 'F-43m (No. 216) holds at 0.2 Å.' }]);
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

    it('bond-angle checks: a window that misses most central atoms warns', () => {
        const { checks } = summarizeBondAngles({
            triplet: ['Se', 'Ta', 'Se'], angleCount: 2, apexCount: 10, coordination: [8, 0, 2],
            binCenters: [89.5, 90.5], counts: [1, 1], meanAngle: 90, stdAngle: 0.5,
            lengths12: { meanLength: 2.6, uniqueBonds: 4 }, sharedEnds: true
        });
        expect(checks).toEqual([{ level: 'warn', text: '80 % of the Ta atoms have no Se–Ta bond in the window: it misses most of the shell.' }]);
    });

    it('symmetry checks: not analysed, undetermined, lower bound', () => {
        expect(summarizeSymmetry({ skipped: true, reason: 'too many', spaceGroup: 'not analysed', orbits: [] }, 0.2).checks[0])
            .toEqual({ level: 'warn', text: 'Not analysed: too many' });
        expect(summarizeSymmetry({ spaceGroup: 'undetermined', orbits: [] }, 0.2).checks[0].level).toBe('note');
        expect(summarizeSymmetry({ spaceGroup: '≥ -43m', orbits: [] }, 0.2).checks[0].text).toMatch(/lower bound/);
        expect(summarizeSymmetry(null, 0.2).checks[0].level).toBe('fail');
    });
});
