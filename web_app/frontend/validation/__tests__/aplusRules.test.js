// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// Arm A+ reads the app's rules (architecture candidate 8, 2026-10-19):
// known-answer tests for each rule as A+ now applies it — the split site
// (isSplitSite, with list_sites' static displacement), the cluster elements
// (the evidence gate's clusterElementsOf) and the lower group's tolerance
// (lowerGroupTolerance, NOISE_FLOOR_MARGIN × the noise floor). Synthetic
// clouds, catalogues and structures, and the DEV cases only: the held-out
// test cases are never built here.

import { beforeAll, describe, expect, it } from 'vitest';
import { buildCase, DEV_CASES, mulberry32 } from '../cases';
import { loadRun } from '../runFolder';
import { appSplitSites, ARM_A_RULES, runArmSweep, SWEEP_RULES } from '../armA';
import {
    directionCandidate, directionResolved, directionTier, isDiffuseOrder, LOW_INDEX_CONE_DEG, neighbourFrequency, splitDisplacement, summarizeSites
} from '../../src/copilotTools';
import { buildClarifier } from '../../src/llm/copilot/clarify';
import { createVerification } from '../../src/llm/copilot/verification';
import { lowerGroupTolerance, NOISE_FLOOR_MARGIN, symmetryNoiseFloor } from '../../src/llm/context/runContext';
import { CLUSTER_CN_MAX, clusterElementsOf, shellCatalogue } from '../../src/llm/context/shells';
import { describeSymmetry } from '../../src/symmetryModel';
import { isSplitSite } from '../../src/workers/clusters';
import { siteEllipsoids } from '../../src/workers/pcaKde';

// --- synthetic site clouds through the engine's site rows ------------------------------

const normal = (rng) => {
    let u = 0;
    while (u === 0) u = rng();
    return Math.sqrt(-2 * Math.log(u)) * Math.cos(2 * Math.PI * rng());
};
// One site's cloud (Cartesian Å, cubic cell, the split or long axis along z):
// equal wells at ±delta (0: none) with Gaussian widths sigma = [x, y, z].
const cloud = ({ n, sigma, delta = 0 }, seed) => {
    const rng = mulberry32(seed);
    return Array.from({ length: n }, () => {
        const well = delta ? (rng() < 0.5 ? delta : -delta) : 0;
        return [sigma[0] * normal(rng), sigma[1] * normal(rng), well + sigma[2] * normal(rng)];
    });
};
const CELL = [[4, 0, 0], [0, 4, 0], [0, 0, 4]];
const rowOf = (displacements) => siteEllipsoids([{
    referenceNumber: 1, element: 'F', elementCounts: { F: displacements.length }, mixed: false,
    count: displacements.length, displacements, siteFractional: [0, 0, 0], copiesPerCell: null
}])[0];
const payloadOf = (rows) => ({ sites: rows, totalAtoms: rows.reduce((sum, row) => sum + row.count, 0), latticeVectors: CELL, supercell: [1, 1, 1], parseWarning: null });
// What A+'s split step needs of a run: the 'sites' engine and the cell.
const fakeRun = (rows) => ({ engine: async () => payloadOf(rows), structure: { latticeVectors: CELL, supercell: [1, 1, 1] } });
// Arm A's frozen rule (ARM_A_RULES), which A+ used until candidate 8.
const armARule = (row) => row.rms[0] / row.rms[1] >= ARM_A_RULES.PROLATE_RATIO && Boolean(row.axisResolved?.[0])
    && Number.isFinite(row.excessKurtosis?.[0]) && row.excessKurtosis[0] <= ARM_A_RULES.SPLIT_KURTOSIS;

describe('A+ reads the app\'s rules (candidate 8): one definition each', () => {
    it('the rules A+ calls are the app\'s own functions', () => {
        expect(SWEEP_RULES.splitSite).toBe(isSplitSite);
        expect(SWEEP_RULES.staticDisplacement).toBe(splitDisplacement);
        expect(SWEEP_RULES.clusterElements).toBe(clusterElementsOf);
        expect(SWEEP_RULES.lowerGroupTolerance).toBe(lowerGroupTolerance);
        // Plan item 2.5: a peak cell resolves a direction within LOW_INDEX_CONE_DEG.
        expect(SWEEP_RULES.directionResolved).toBe(directionResolved);
        // The displacement-anisotropy interview (2026-10-24): the tool's direction tiers.
        expect(SWEEP_RULES.directionTier).toEqual({ tier: directionTier, candidate: directionCandidate, neighbour: neighbourFrequency });
        // The local-order interview (2026-10-24): the diffuse verdict is order in antiphase only, the tool's facts' reading.
        expect(SWEEP_RULES.diffuseOrder).toBe(isDiffuseOrder);
        expect(NOISE_FLOOR_MARGIN).toBe(1.25);
        expect(CLUSTER_CN_MAX).toBe(4.5);
        expect(LOW_INDEX_CONE_DEG).toBe(15);
        expect([14.9, 15, 15.1, 37.4, null].map(directionResolved)).toEqual([true, true, false, false, false]);
    });
});

describe('the split site as A+ applies it (isSplitSite, list_sites\' δ)', () => {
    it('a prolate Gaussian is never split, where arm A\'s rule called 3–14 % of them split at 100–216 copies', async () => {
        const SEEDS = 1000;
        for (const [n, armAAtLeast] of [[216, 0.02], [100, 0.08]]) {
            let armA = 0;
            for (let seed = 1; seed <= SEEDS; seed += 1) {
                const row = rowOf(cloud({ n, sigma: [0.05, 0.05, 0.09] }, 104729 * seed + n));
                if (armARule(row)) armA += 1;
                const { found } = await appSplitSites(fakeRun([row]), []);
                expect(found.size, `n ${n} seed ${seed}`).toBe(0);
            }
            // Measured 3.6 % at 216 copies and 14 % at 100 (2,000 seeds each,
            // 2026-10-19): the false calls the app's 3-error bar removes.
            expect(armA / SEEDS).toBeGreaterThan(armAAtLeast);
        }
    }, 120000);

    it('a planted two-site split (±0.15 Å over a 0.05 Å thermal width, 216 copies) is found, with list_sites\' δ, which reads back the planted δ', async () => {
        const SEEDS = 300;
        let covered = 0;
        for (let seed = 1; seed <= SEEDS; seed += 1) {
            const row = rowOf(cloud({ n: 216, sigma: [0.05, 0.05, 0.05], delta: 0.15 }, 7919 * seed + 3));
            const readings = [];
            const { found } = await appSplitSites(fakeRun([row]), readings);
            expect(found.has('F'), `seed ${seed}`).toBe(true);
            const { magnitude } = found.get('F');
            // The δ list_sites gives the model, to its rounding.
            const [listed] = summarizeSites(payloadOf([row])).data.sites;
            expect(Number(magnitude.toPrecision(3))).toBe(listed.delta_A);
            expect(readings[0]).toMatchObject({ page: 'list_sites', value: { delta_A: listed.delta_A, delta_se_A: listed.delta_se_A } });
            const { se } = splitDisplacement(row);
            if (Math.abs(magnitude - 0.15) <= 1.96 * se) covered += 1;
        }
        expect(covered / SEEDS).toBeGreaterThan(0.9);
        expect(covered / SEEDS).toBeLessThan(0.99);
    }, 120000);

    it('one claim per element: the split site with the largest δ', async () => {
        const rows = [0.1, 0.2, 0.15].map((delta, index) => ({
            ...rowOf(cloud({ n: 512, sigma: [0.05, 0.05, 0.05], delta }, 11 + index)), referenceNumber: index + 1
        }));
        const { found } = await appSplitSites(fakeRun(rows), []);
        expect([...found.keys()]).toEqual(['F']);
        expect(found.get('F').site.referenceNumber).toBe(2);
    });
});

// The elements the evidence gate (the verification) asks analyze_clusters for
// before "no order", read off its finding with nothing run.
const gateClusterElements = (catalog) => createVerification({ catalog }).judge('There is no local order beyond thermal motion.').done.audit[0].missing
    .filter((item) => item.startsWith('analyze_clusters for ')).map((item) => item.split(' ')[2]);

describe('the cluster elements as A+ applies them (the evidence gate\'s rule)', () => {
    it('an element whose own first shell holds 0.5–4.5 neighbours: dimers, trimers, tetrahedra and octahedra, not networks', () => {
        const cns = { D: 1, T: 2, Te: 3, O: 4, B: 4.5, Ic: 5, S: 6, P: 8, F: 12, I: 0.3, U: null };
        const catalog = { elements: Object.keys(cns).sort(), shell: (a, b) => (a === b ? { cn: cns[a] } : null) };
        expect(clusterElementsOf(catalog)).toEqual(['B', 'D', 'O', 'T', 'Te']);
        expect(gateClusterElements(catalog)).toEqual(clusterElementsOf(catalog));
        expect(clusterElementsOf(null)).toEqual([]);
    });

    describe('on the dev cases (planted dimers, trimers and octahedra)', () => {
        const runs = new Map();
        beforeAll(async () => {
            for (const item of DEV_CASES) runs.set(item.id, await loadRun(buildCase(item), { ladder: false }));
        }, 120000);

        it('A+ tries clusters for the gate\'s elements: the planted V₂ dimers, Mo₃ triangles, W₆ octahedra and the F around Cu and Ag, and nothing else', async () => {
            const tried = {};
            for (const [id, run] of runs) {
                const catalogue = shellCatalogue(run.props);
                const { readings } = await runArmSweep(run);
                const elements = readings.filter((reading) => reading.page === 'analyze_clusters').map((reading) => reading.what.split(' ')[0]);
                expect(elements, id).toEqual(gateClusterElements(catalogue));
                if (elements.length) tried[id] = elements;
            }
            expect(tried).toEqual({
                'dev-breathing': ['Mo'], 'dev-jahn-teller': ['F'], 'dev-cell-wrap': ['V'], 'dev-jt-off-centre': ['F'], 'dev-octahedral-clusters': ['W']
            });
        }, 120000);

        it('the elements the old 3.6 Å rule also tried are extended networks the engine refuses (CN 6–12)', async () => {
            let refused = 0;
            for (const [id, run] of runs) {
                const catalogue = shellCatalogue(run.props);
                const kept = new Set(clusterElementsOf(catalogue));
                for (const shell of catalogue?.pairs || []) {
                    if (shell.centre !== shell.partner || !shell.window || shell.peaks[0].r > 3.6 || kept.has(shell.centre)) continue;
                    expect(shell.cn, `${id} ${shell.centre}`).toBeGreaterThan(CLUSTER_CN_MAX);
                    const result = await run.engine('clusters', { element: shell.centre, rMin: shell.window[0], rMax: shell.window[1] });
                    expect(result.verdict.window.extendedNetwork.because, `${id} ${shell.centre}`).not.toBe('finite_clusters');
                    refused += 1;
                }
            }
            // O in seven perovskites, Fe and Si in B20 FeSi, Mg and O in rock salt.
            expect(refused).toBe(11);
        }, 120000);
    });
});

describe('the lower group\'s tolerance as A+ applies it (NOISE_FLOOR_MARGIN × the noise floor)', () => {
    // A KTaO₃-like cell with Ta moved δ along c, and site means whose noise
    // floor is 3√2 · 0.1/√450 = 0.020 Å: P4mm holds to about 1.6 δ.
    const perovskite = (delta) => ({
        latticeVectors: CELL, supercell: [1, 1, 1], elementCounts: { K: 1, Ta: 1, O: 3 },
        basis: [
            ['K', [0, 0, 0]], ['Ta', [0.5, 0.5, 0.5 + delta / 4]], ['O', [0.5, 0.5, 0]], ['O', [0.5, 0, 0.5]], ['O', [0, 0.5, 0.5]]
        ].map(([el, frac], index) => ({ el, frac, referenceNumber: index + 1, dispA: 0.1, count: 450 }))
    });

    it('is 1.25 × the floor at two digits, and 0.025 Å without counts', () => {
        expect(symmetryNoiseFloor(perovskite(0).basis)).toBeCloseTo(0.02, 6);
        expect(lowerGroupTolerance(perovskite(0).basis)).toBe(0.025);
        expect(lowerGroupTolerance([{ dispA: 0.2, count: 400 }])).toBe(0.053);
        expect(lowerGroupTolerance([{ el: 'K', frac: [0, 0, 0] }])).toBe(0.025);
    });

    it('a lowering whose group holds only to 1.2 × the floor is not read; one that holds to 2 × is, as the Workbench offers it', () => {
        const floor = 0.02;
        // δ = 0.0147 Å: P4mm holds to 0.0235 Å (1.18 × the floor).
        const marginal = perovskite(0.0147);
        expect(describeSymmetry(marginal, floor).spaceGroup).toBe('P4mm');
        expect(describeSymmetry(marginal, SWEEP_RULES.lowerGroupTolerance(marginal.basis)).spaceGroup).toBe('Pm-3m');
        // δ = 0.027 Å: P4mm holds to 0.043 Å (2.2 × the floor).
        const clear = perovskite(0.027);
        expect(describeSymmetry(clear, SWEEP_RULES.lowerGroupTolerance(clear.basis)).spaceGroup).toBe('P4mm');
        // The Workbench's question counts the same groups.
        const asks = (to) => buildClarifier({ symmetry: { ladder: [{ spaceGroup: 'P4mm', from: 0, to }, { spaceGroup: 'Pm-3m', from: to, to: 1 }] } }, { noiseFloor: floor })
            .beforeTurn('What is the space group?');
        expect(asks(0.0235)).toBeNull();
        expect(asks(0.043)?.request.options.map((option) => option.label)).toEqual(['Pm-3m', 'P4mm']);
    });

    it('on the dev cases A+ reads the truth\'s group at the tool\'s default tolerance, every lower rung below the floor', async () => {
        for (const item of DEV_CASES) {
            const built = buildCase(item);
            const run = await loadRun(built, { ladder: true });
            const { basis } = run.props.structure;
            const floor = symmetryNoiseFloor(basis);
            const tolerance = lowerGroupTolerance(basis);
            const lowering = run.tools.find((tool) => tool.name === 'measure_symmetry_lowering');
            expect(lowering.parameters.properties.tolerance_A.description, item.id).toBe(`Tolerance of the lower group (Å). Default ${tolerance}.`);
            const { readings, verdict } = await runArmSweep(run);
            expect(readings[0], item.id).toMatchObject({ what: `space group at ${tolerance} Å`, noise_floor_A: Number(floor.toPrecision(3)) });
            // The truth's group (dev-breathing's truth scores no symmetry: its
            // configuration's own group, R3m).
            const group = built.truth.symmetry?.space_group ?? 'R3m';
            expect(verdict.space_group, item.id).toBe(group);
            // The noise rungs end well below the floor (at most 0.67 × it on
            // dev), so the margin moves no dev group.
            const ladder = run.props.symmetry.ladder;
            const truth = ladder.findIndex((brick) => brick.spaceGroup === group);
            expect(truth, item.id).toBeGreaterThanOrEqual(0);
            ladder.slice(0, truth).forEach((brick) => expect(brick.to, `${item.id} ${brick.spaceGroup}`).toBeLessThan(0.7 * floor));
            expect(ladder[truth].to, item.id).toBeGreaterThan(tolerance);
        }
    }, 120000);
});
