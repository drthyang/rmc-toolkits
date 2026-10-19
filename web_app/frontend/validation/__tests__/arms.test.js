// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// Arm A (the default pages) and A+ (the engine sweep) on the DEV cases:
// pinned, because they are the bar a Copilot run must beat, and a change to
// an engine or a page default that moves them must be seen. The test set's
// arm-A numbers come from the CLI (`arm-a --split test --final`), never here.

import { describe, expect, it } from 'vitest';
import { shellCatalogue } from '../../src/llm/context/shells';
import { buildCase, DEV_CASES } from '../cases';
import { loadRun } from '../runFolder';
import { ARM_A_RULES, runArmA, runArmSweep } from '../armA';
import { scoreVerdict } from '../scoring';

const runBoth = async (item) => {
    const built = buildCase(item);
    const run = loadRun(built, { ladder: false });
    const a = await runArmA(run);
    const sweep = await runArmSweep(run);
    return { built, run, a, sweep, scoreA: scoreVerdict(built.truth, a.verdict), scoreSweep: scoreVerdict(built.truth, sweep.verdict) };
};

// [case, arm A score, arm A spurious claims, A+ score, A+ spurious claims]
// Batch 1 (2026-10-17): the sweep gained measure_symmetry_lowering and the
// general-wavevector test (review L4), and dev-symmetry-lowering's A+ went
// from 0.8333 to 1 (the shift is named); every other pin is unchanged.
// Plan item 4.6: dev-fit-unconverged (dev only), scored on its
// structure (no order, Fm-3m); its fit files change nothing the arms read.
// Plan item 1.6 ii: the judgment subset (dev only), the last five; every pin
// above them is unchanged.
// Candidate 8 (2026-10-19, A+ reads the app's rules): no pin moved; only A+'s
// readings (and dev-offcentre's δ, 0.201 → 0.202 Å) did, pinned below.
// Plan item 2.5 (protocol 2026-10-19.16, the Displacement Directions page at
// the engine's resolution, its readouts the raw histogram's): one pin moved,
// dev-jt-off-centre's A and A+ 0.8571 → 0.7857. At its 216 copies Auto is
// ν = 1, twelve cells with [0 0 1] on a boundary, so the now significant Ag
// peak names [0 −0.62 1], 31.7° off the planted axis (pinned below).
const BAR = [
    ['dev-offcentre', 1, [], 1, []],
    ['dev-zone-boundary', 0.6667, ['off_centre Ti'], 1, []],
    ['dev-breathing', 0.5, [], 1, []],
    ['dev-jahn-teller', 0.9167, [], 0.9167, []],
    ['dev-cell-wrap', 0.8333, [], 1, []],
    ['dev-thermal-null', 1, [], 1, []],
    ['dev-symmetry-lowering', 0.3333, [], 1, []],
    ['dev-chirality-trap', 1, [], 1, []],
    ['dev-chiral', 1, [], 1, []],
    ['dev-fit-unconverged', 1, [], 1, []],
    ['dev-jt-off-centre', 0.7857, [], 0.7857, []],
    ['dev-eight-site', 0.6667, ['jahn_teller Nb'], 0.6667, ['jahn_teller Nb']],
    ['dev-octahedral-clusters', 0.8333, [], 1, []],
    ['dev-weak-order-a', 0.8333, [], 1, []],
    ['dev-weak-order-b', 0.8333, [], 0.8333, []]
];

describe('the model-free arms on the dev cases (the bar)', () => {
    it.each(BAR)('%s: arm A %s, A+ %s', async (id, scoreA, spuriousA, scoreSweep, spuriousSweep) => {
        const { scoreA: a, scoreSweep: s } = await runBoth(DEV_CASES.find((item) => item.id === id));
        expect(a.score).toBe(scoreA);
        expect(a.spurious.map((c) => `${c.kind} ${c.element}`)).toEqual(spuriousA);
        expect(s.score).toBe(scoreSweep);
        expect(s.spurious.map((c) => `${c.kind} ${c.element}`)).toEqual(spuriousSweep);
    });

    it('where the default pages mislead, the readings say which page and value', async () => {
        const lowering = await runBoth(DEV_CASES.find((item) => item.id === 'dev-symmetry-lowering'));
        // The Detected SG card at 0.2 Å keeps the parent group; the sweep reads at the noise floor.
        expect(lowering.a.verdict).toMatchObject({ space_group: 'Pm-3m', centrosymmetric: true, polar: false });
        expect(lowering.a.readings[0]).toMatchObject({ page: 'Detected SG card', value: 'Pm-3m' });
        expect(lowering.sweep.verdict).toMatchObject({ space_group: 'P4mm', centrosymmetric: false, polar: true });
        // The shift tool names the shift the card hides (batch 1, review L4).
        const [shift] = lowering.sweep.verdict.features;
        expect(shift).toMatchObject({ kind: 'symmetry_lowering', element: 'Ta', direction_uvw: [0, 0, 1] });
        expect(Math.abs(shift.magnitude_A - 0.07)).toBeLessThan(0.01);

        const zone = await runBoth(DEV_CASES.find((item) => item.id === 'dev-zone-boundary'));
        const ti = zone.a.readings.find((r) => r.page === 'PCA Ellipsoid' && r.what.endsWith('Ti'));
        expect(ti.reading).toMatch(/^split site along PC1/);
        expect(zone.sweep.verdict.features).toEqual([expect.objectContaining({ kind: 'modulation', element: 'Ti', q: [0, 0, 0.5] })]);

        const jt = await runBoth(DEV_CASES.find((item) => item.id === 'dev-jahn-teller'));
        const shell = jt.a.readings.find((r) => r.page === 'Dashboard partial g(r)');
        expect(shell.value.peaks_A).toHaveLength(2);
        expect(jt.a.verdict.features).toContainEqual(expect.objectContaining({ kind: 'jahn_teller', element: 'Cu' }));
    });

    it('A+ reads the cluster and modulation engines\' own verdicts (candidate 5, 2026-10-19)', async () => {
        const breathing = await runBoth(DEV_CASES.find((item) => item.id === 'dev-breathing'));
        const clusters = breathing.sweep.readings.find((r) => r.page === 'analyze_clusters');
        expect(clusters.value).toMatchObject({ breathing: 'true (above_null)', distorted: 'false (within_null)', rigid: 'false (within_null)' });
        expect(breathing.sweep.verdict.features).toEqual([{ kind: 'cluster_breathing', element: 'Mo' }]);
        const zone = await runBoth(DEV_CASES.find((item) => item.id === 'dev-zone-boundary'));
        const modulations = zone.sweep.readings.find((r) => r.page === 'analyze_modulations');
        expect(modulations.reading).toBe('cell doubling: Ti at (0, 0, 0.5); order at a general q: none; long wave: Sr at (0.125, 0, 0)');
        const jt = await runBoth(DEV_CASES.find((item) => item.id === 'dev-jahn-teller'));
        expect(jt.sweep.readings.find((r) => r.page === 'analyze_clusters').value)
            .toMatchObject({ distorted: 'true (axis_family)', splitMemberSites: 'true (half_or_more_split)', rigid: 'false (split_sites)' });
    });

    it('A+ reads a modulation\'s direction from the engine\'s polarisation, not the PC1 of a prolate site (plan item 4.3)', async () => {
        // dev-zone-boundary: the Ti order's polarisation (5.5° from the planted
        // [0 0 1], 5.9° as the claim's two decimals score; the most prolate
        // site's PC1 it read before, [−0.03 0.1 1], scored 6.0°).
        const zone = await runBoth(DEV_CASES.find((item) => item.id === 'dev-zone-boundary'));
        expect(zone.sweep.verdict.features).toEqual([{ kind: 'modulation', element: 'Ti', q: [0, 0, 0.5], direction_uvw: [-0.05, 0.09, 1] }]);
        expect(zone.scoreSweep.features[0]).toMatchObject({ direction_ok: true, q_ok: true, accurate: true });
        expect(zone.scoreSweep.features[0].direction_error_deg).toBeCloseTo(5.88, 2);
        const reading = zone.sweep.readings.find((r) => r.page === 'analyze_modulations');
        expect(reading.value.Ti.polarisation).toEqual([{ q: '(0, 0, 0.5)', uvw: [-0.05, 0.09, 1], share: 1 }]);
        expect(reading.value.Sr.polarisation).toBeUndefined();
        // dev-weak-order-a: no Nb site is prolate (the order is far below the
        // thermal width), so the claim had no direction; the polarisation gives
        // one (10° from the planted axis, which the truth does not score).
        const weak = await runBoth(DEV_CASES.find((item) => item.id === 'dev-weak-order-a'));
        expect(weak.sweep.verdict.features).toEqual([{ kind: 'modulation', element: 'Nb', q: [0.5, 0.5, 0], direction_uvw: [0.17, 0.06, 1] }]);
    });

    it('A+ reads the app\'s rules: the split site, the cluster elements, the lower group\'s tolerance (candidate 8, 2026-10-19)', async () => {
        // The split site: isSplitSite, with list_sites' δ. On dev-jahn-teller
        // two F sites (κ₁ −0.95 and −0.93 at 216 copies: 2.8 errors) are split
        // by arm A's rule, not by the app's; the F claim (the largest δ) stays.
        const jt = await runBoth(DEV_CASES.find((item) => item.id === 'dev-jahn-teller'));
        const site = (arm, n) => arm.readings.find((r) => /^(PCA Ellipsoid|list_sites)$/.test(r.page) && r.what === `site ${n} F`);
        for (const n of [5, 6]) {
            expect(site(jt.a, n).reading).toMatch(/^split site along PC1/);
            expect(site(jt.sweep, n)).toMatchObject({ page: 'list_sites', reading: 'not a split site' });
            expect(site(jt.sweep, n).value.pc1_kurtosis_z).toBeGreaterThan(-3);
        }
        expect(site(jt.sweep, 3).value).toMatchObject({ delta_A: 0.16, delta_se_A: 0.0043, pc1_kurtosis_z: -4.43 });
        // Its direction since plan item 2.5: the map's peak on the engine's
        // ν = 1 tiling (significant there; on the old default it was not, and
        // the claim took PC1, [1 0 −0.01]). F off-centring is tolerated, unscored.
        expect(jt.sweep.verdict.features[0]).toEqual({ kind: 'off_centre', element: 'F', direction_uvw: [-1, 0, -0.62], magnitude_A: 0.16 });
        const offcentre = await runBoth(DEV_CASES.find((item) => item.id === 'dev-offcentre'));
        // ν = 2 holds [0 0 1] as a cell centre (on ν = 10 the peak read [0.1 0.06 −1]).
        expect(offcentre.sweep.verdict.features).toEqual([{ kind: 'off_centre', element: 'Nb', direction_uvw: [0, 0, -1], magnitude_A: 0.202 }]);
        expect(offcentre.a.verdict.features[0].magnitude_A).toBe(0.201);
        // The cluster elements: the gate's (CN 0.5–4.5); perovskite O (CN 8) no longer tried.
        expect(offcentre.sweep.readings.filter((r) => r.page === 'analyze_clusters')).toEqual([]);
        // The lower group: 1.25 × the noise floor (0.026 Å → 0.032 Å), the tool's default.
        const lowering = await runBoth(DEV_CASES.find((item) => item.id === 'dev-symmetry-lowering'));
        expect(lowering.sweep.readings[0]).toMatchObject({ page: 'detect_symmetry at 1.25 × the noise floor', what: 'space group at 0.032 Å', value: 'P4mm', noise_floor_A: 0.026 });
        expect(lowering.sweep.readings.find((r) => r.page === 'measure_symmetry_lowering').reading).toContain('(the lower group P4mm at 0.032 Å)');
    });

    // The room a model has (plan 1.6 ii): on each judgment case, which of
    // A+'s fixed rules misses or misreads the planted truth, read from its
    // own readings, and that the engines hold the right answer.
    it('the judgment subset: where A+\'s fixed joining rules fail, the readings show which rule', async () => {
        const judged = (id) => runBoth(DEV_CASES.find((item) => item.id === id));
        const found = (score) => score.features.map((feature) => `${feature.planted ?? feature.kind}:${feature.found}`);

        // Two features: Ag is split, and by rule a split centre explains its
        // own split shell, so the Jahn–Teller elongation is never claimed.
        const two = await judged('dev-jt-off-centre');
        expect(two.sweep.readings.find((r) => r.page === 'Dashboard partial g(r)').value.peaks_A).toHaveLength(2);
        expect(two.sweep.verdict.features.some((feature) => feature.kind === 'jahn_teller')).toBe(false);
        expect(found(two.scoreSweep)).toEqual(['jahn_teller:false', 'off_centre:true']);
        // A side miss since plan item 2.5: the Ag direction map at the engine's
        // resolution, ν = 1 for 216 copies, whose twelve cell centres sit 31.7°
        // from every cube axis. Its peak is significant, so A and A+ take its
        // direction ([0 −0.62 1]) over PC1's (within 2° of [0 0 1]).
        const agMap = two.sweep.readings.find((r) => r.page === 'Displacement Directions' && r.what.endsWith('Ag'));
        expect(agMap.value).toMatchObject({ frequency: 1, peak_uvw: [0, -0.62, 1] });
        expect(agMap.value.peak_sigma).toBeGreaterThanOrEqual(ARM_A_RULES.SIGMA);
        const ag = two.scoreSweep.features.find((feature) => feature.planted === 'off_centre');
        expect(ag).toMatchObject({ direction_ok: false, magnitude_ok: true });
        expect(ag.direction_error_deg).toBeCloseTo(31.7, 0);

        // Two readings: the eight wells average to a round Nb ellipsoid, so
        // no split site; the 3 + 3 shell it leaves unexplained becomes a
        // Jahn–Teller claim. The direction map at a finer tiling has the
        // <111> lobe that the default one does not resolve.
        // (Since candidate 8 A+ reads list_sites' split rule: the round
        // ellipsoid resolves no PC1, so no split site, as arm A's page reads it.)
        const eight = await judged('dev-eight-site');
        const nb = eight.sweep.readings.find((r) => r.page === 'list_sites' && r.what.endsWith('Nb'));
        expect(nb.reading).toBe('not a split site');
        expect(nb.value.pc1_resolved).toBe(false);
        const nbPage = eight.a.readings.find((r) => r.page === 'PCA Ellipsoid' && r.what.endsWith('Nb'));
        expect(nbPage.reading).toBe('thermal ellipsoid');
        expect(nbPage.value.pc1_over_pc2).toBeLessThan(ARM_A_RULES.PROLATE_RATIO);
        expect(found(eight.scoreSweep)).toEqual(['off_centre:false']);
        const engine = (await eight.run).engine;
        const coarse = await engine('orientation', { referenceNumber: 2 });
        const fine = await engine('orientation', { referenceNumber: 2, frequency: 3 });
        expect(coarse.peakSignificance).toBeLessThan(ARM_A_RULES.SIGMA);
        expect(coarse.mapSignificance).toBeGreaterThan(ARM_A_RULES.SIGMA);
        expect(fine.peakSignificance).toBeGreaterThan(ARM_A_RULES.SIGMA);
        expect(fine.peakDirection.map((value) => Math.abs(value))).toEqual(fine.peakDirection.map(() => expect.closeTo(1 / Math.sqrt(3), 2)));

        // A choice only the user can make, the one control: since the one
        // shell reading (protocol 2026-10-19.8) A+ reads the Workbench's own
        // joined W–W shell (CN 4, its default) and finds the distorted
        // octahedra. The first peak alone, the window the Workbench asks
        // about when a call cuts the shell, holds no cluster. The judgment
        // (stating the default and what the other window gives) shows in the
        // answer, not in the score.
        const octahedra = await judged('dev-octahedral-clusters');
        expect(found(octahedra.scoreSweep)).toEqual(['cluster_distortion:true']);
        const run = await octahedra.run;
        const shell = shellCatalogue(run.props).shell('W', 'W');
        expect(shell.joined).toBe(true);
        expect(shell.firstPeak.changesCn).toBe(true);
        expect(octahedra.sweep.readings.find((r) => r.page === 'analyze_clusters').what)
            .toBe(`W clusters in ${shell.window[0]}–${shell.window[1]} Å (engine verdict)`);
        const cut = await run.engine('clusters', { element: 'W', rMin: shell.firstPeak.window[0], rMax: shell.firstPeak.window[1] });
        expect(cut.verdict.window.extendedNetwork.because).toBe('no_clusters');

        // A signal at the edge of significance: one weak order, a seed on
        // each side of the fixed 3σ bar; A+ claims the one above only.
        for (const [id, ordered] of [['dev-weak-order-a', true], ['dev-weak-order-b', false]]) {
            const weak = await judged(id);
            const payload = await (await weak.run).engine('modulations', {});
            const point = payload.elements.find((entry) => entry.element === 'Nb').zoneBoundary.find((entry) => entry.q.join() === '0.5,0.5,0');
            expect(point.sigma, id).toBeGreaterThanOrEqual(2.5);
            expect(point.sigma, id).toBeLessThanOrEqual(4);
            expect(point.ordered, id).toBe(ordered);
            expect(found(weak.scoreSweep), id).toEqual([`modulation:${ordered}`]);
        }
    }, 120000);
});
