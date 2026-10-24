// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// A direction's tier (the displacement-anisotropy interview, 2026-10-24): a
// finding at ≥ 3σ (unrounded) and ≥ 1.2× isotropic that holds at the
// neighbouring resolution; a hint at ≥ 3σ under 1.2×, 2.5–3σ, or held only at
// the engine's own resolution; nothing below 2.5σ. Known answers on the rule's
// edges, and on the Demo and a dev case through the tool, as the app runs it.

import { describe, expect, it } from 'vitest';
import {
    buildCopilotTools, DIRECTION_FINDING_ENHANCEMENT, DIRECTION_HINT_SIGMA, directionCandidate, directionHolds, directionTier, neighbourFrequency,
    SIGNIFICANCE_SIGMA
} from '../copilotTools';
import { buildCase, DEV_CASES } from '../../validation/cases';
import { loadRun } from '../../validation/runFolder';
import { loadDemoAssistantRun } from './modelVisible/helpers';

// An engine payload's peak: its σ, enhancement, direction and cell radius.
const peak = (peakSignificance, peakEnhancement = 2, peakDirection = [0, 0, 1], peakCellRadiusDeg = 20.9) => ({
    peakSignificance, peakEnhancement, peakDirection, peakCellRadiusDeg
});

describe('the tiers on their edges', () => {
    it('the constants are the user\'s', () => {
        expect([SIGNIFICANCE_SIGMA, DIRECTION_FINDING_ENHANCEMENT, DIRECTION_HINT_SIGMA]).toEqual([3, 1.2, 2.5]);
        expect([1, 2, 5, 24].map(neighbourFrequency)).toEqual([2, 1, 4, 23]);
    });

    it('nothing below 2.5σ, a hint to 3σ, and at 3σ a hint under 1.2× isotropic', () => {
        expect(directionTier(peak(2.4999))).toEqual({ tier: 'none', because: 'below_hint_sigma' });
        expect(directionTier(peak(null))).toEqual({ tier: 'none', because: 'below_hint_sigma' });
        expect(directionTier(peak(2.5))).toEqual({ tier: 'hint', because: 'below_sigma' });
        // Unrounded: 2.996σ is no 3σ.
        expect(directionTier(peak(2.996))).toEqual({ tier: 'hint', because: 'below_sigma' });
        expect(directionTier(peak(3, 1.1999), peak(5))).toEqual({ tier: 'hint', because: 'small_enhancement' });
        expect(directionCandidate(peak(3, 1.1999))).toBe(false);
        expect(directionCandidate(peak(3, 1.2))).toBe(true);
    });

    it('a candidate is a finding only where it holds at the neighbouring resolution: ≥ 3σ there, the same axis within the two cells\' radii', () => {
        const first = peak(4, 2, [0, 0, 1], 20.9);
        expect(directionTier(first, peak(3.2, 1.5, [Math.sin(0.5), 0, Math.cos(0.5)], 37.4))).toEqual({ tier: 'finding', because: 'holds' });
        // The opposite well of a split site: the same axis.
        expect(directionTier(first, peak(9, 3, [0, 0, -1], 37.4))).toEqual({ tier: 'finding', because: 'holds' });
        expect(directionHolds(first, peak(9, 3, [0, 0, -1], 37.4))).toEqual({ holds: true, apartDeg: 0 });
        // Below 3σ there, or another axis: held only at the engine's resolution.
        expect(directionTier(first, peak(2.9, 1.5, [0, 0, 1], 37.4))).toEqual({ tier: 'hint', because: 'not_held' });
        expect(directionTier(first, peak(6, 2, [1, 0, 0], 37.4))).toEqual({ tier: 'hint', because: 'not_held' });
        expect(directionHolds(first, peak(6, 2, [1, 0, 0], 37.4)).apartDeg).toBeCloseTo(90, 9);
        // Unchecked, or read at another resolution than the engine's: a hint at most.
        expect(directionTier(first)).toEqual({ tier: 'hint', because: 'unchecked' });
        expect(directionTier(first, peak(9), { auto: false })).toEqual({ tier: 'hint', because: 'not_auto' });
    });
});

describe('the tiers through the tool (analyze_displacement_directions)', () => {
    // The run's engine, counting its orientation requests.
    const counted = (engine) => {
        const asked = [];
        const request = (kind, params, options) => {
            if (kind === 'orientation') asked.push(params.frequency ?? 'auto');
            return engine(kind, params, options);
        };
        return { request, asked };
    };

    it('the Demo: Ta site 17 is a finding, Ta site 25 a hint it does not hold, site 1 a hint below 3σ, Se site 2 nothing; the second call only for a candidate', async () => {
        const { copilot, run } = await loadDemoAssistantRun();
        const { request, asked } = counted(copilot.engine);
        const tool = buildCopilotTools({ request, structure: run.structure, elements: Object.keys(run.structure.elementCounts).sort() })
            .find((item) => item.name === 'analyze_displacement_directions');
        const read = async (ref) => {
            asked.length = 0;
            const result = await tool.run({ site_ref: ref });
            return { ...result, asked: [...asked] };
        };
        const finding = await read(17);
        expect(finding.data.peak).toMatchObject({ tier: 'finding', significance_sigma: 3.15 });
        expect(finding.data.neighbour).toEqual({ frequency: 1, significance_sigma: 3.83, enhancement: 1.52, radius_deg: 37.4, apart_deg: 31.7, holds: true });
        expect(finding.asked).toEqual(['auto', 1]);
        expect(finding.checks.find((check) => check.level === 'ok').text).toMatch(/^Preferred direction \(a finding\): the peak cell is 1\.95× isotropic at 3\.15σ/);
        const held = await read(25);
        expect(held.data.peak).toMatchObject({ tier: 'hint', significance_sigma: 4.47 });
        expect(held.data.neighbour).toMatchObject({ frequency: 1, significance_sigma: 2.33, holds: false });
        expect(held.facts.direction).toMatchObject({ tier: 'hint', because: 'not_held' });
        const below = await read(1);
        expect(below.data.peak).toMatchObject({ tier: 'hint', significance_sigma: 2.51 });
        expect(below.data).not.toHaveProperty('neighbour');
        expect(below.asked).toEqual(['auto']);
        expect(below.checks.some((check) => /^A hint of a direction, not a finding: the peak is 2\.51σ, between 2\.5σ and 3σ/.test(check.text))).toBe(true);
        const none = await read(2);
        expect(none.data.peak.tier).toBe('none');
        expect(none.asked).toEqual(['auto']);
        // Beside the engine's resolution (Auto − 1, the call the skill asks for): the other of the pair is the Auto ν itself.
        asked.length = 0;
        const beside = await tool.run({ site_ref: 17, frequency: 1 });
        expect(asked).toEqual([1, 2]);
        expect(beside.data.neighbour).toMatchObject({ frequency: 2, significance_sigma: 3.15, holds: true });
        expect(beside.facts.direction).toMatchObject({ tier: 'finding', because: 'holds', frequency: 1 });
        // Two away from it: no second call, a hint at most.
        asked.length = 0;
        const far = await tool.run({ site_ref: 17, frequency: 4 });
        expect(asked).toEqual([4]);
        expect(far.facts.direction.tier === 'finding').toBe(false);
    }, 120000);

    it('a dev split site: dev-jahn-teller\'s F holds on the opposite well at the neighbouring resolution, the same axis', async () => {
        const copilotRun = await loadRun(buildCase(DEV_CASES.find((item) => item.id === 'dev-jahn-teller')));
        const tool = copilotRun.tools.find((item) => item.name === 'analyze_displacement_directions');
        const { data } = await tool.run({ site_ref: 3 });
        expect(data.frequency).toBe(1);
        expect(data.neighbour).toMatchObject({ frequency: 2, holds: true, apart_deg: 31.7 });
        expect(data.peak.tier).toBe('finding');
    }, 120000);
});
