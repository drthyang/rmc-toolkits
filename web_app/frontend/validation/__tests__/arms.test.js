// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// Arm A (the default pages) and A+ (the engine sweep) on the DEV cases:
// pinned, because they are the bar a Copilot run must beat, and a change to
// an engine or a page default that moves them must be seen. The test set's
// arm-A numbers come from the CLI (`arm-a --split test --final`), never here.

import { describe, expect, it } from 'vitest';
import { buildCase, DEV_CASES } from '../cases';
import { loadRun } from '../runFolder';
import { runArmA, runArmSweep } from '../armA';
import { scoreVerdict } from '../scoring';

const runBoth = async (item) => {
    const built = buildCase(item);
    const run = loadRun(built, { ladder: false });
    const a = await runArmA(run);
    const sweep = await runArmSweep(run);
    return { built, a, sweep, scoreA: scoreVerdict(built.truth, a.verdict), scoreSweep: scoreVerdict(built.truth, sweep.verdict) };
};

// [case, arm A score, arm A spurious claims, A+ score, A+ spurious claims]
// Batch 1 (2026-10-17): the sweep gained measure_symmetry_lowering and the
// general-wavevector test (review L4), and dev-symmetry-lowering's A+ went
// from 0.8333 to 1 (the shift is named); every other pin is unchanged.
const BAR = [
    ['dev-offcentre', 1, [], 1, []],
    ['dev-zone-boundary', 0.6667, ['off_centre Ti'], 1, []],
    ['dev-breathing', 0.5, [], 1, []],
    ['dev-jahn-teller', 0.9167, [], 0.9167, []],
    ['dev-cell-wrap', 0.8333, [], 1, []],
    ['dev-thermal-null', 1, [], 1, []],
    ['dev-symmetry-lowering', 0.3333, [], 1, []],
    ['dev-chirality-trap', 1, [], 1, []],
    ['dev-chiral', 1, [], 1, []]
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
});
