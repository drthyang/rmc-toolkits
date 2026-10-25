// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// Plan item 2.5: the Displacement Directions page's settings. The page opens
// at ν = 10 with the map smoothed 4× (the user's choice, 2026-10-24), its
// readouts never carry the smoothing (a display setting), and arm A keeps the
// settings frozen at protocol 2026-10-24.6 (Auto, unsmoothed).

import { describe, expect, it } from 'vitest';
import {
    autoResolutionLabel,
    cellCountOf,
    DIRECTION_DEFAULTS,
    directionMapRequest,
    directionReadoutRequest,
    FREQUENCY_OPTIONS
} from '../directionSettings';
import { ARM_A_RULES } from '../../validation/armA';
import { buildCopilotTools } from '../copilotTools';

describe('the Displacement Directions page settings (plan item 2.5)', () => {
    it('open at ν = 10 with the map smoothed 4×, the readouts at smoothing 0 (the user\'s choice, 2026-10-24)', () => {
        // ν = 10 is a resolution choice (1002 cells); 4× is four neighbour passes.
        expect(DIRECTION_DEFAULTS).toEqual({ frequency: 10, weight: 'count', frame: 'cartesian', minQuantile: 0, smoothing: 4 });
        expect(FREQUENCY_OPTIONS).toContain(DIRECTION_DEFAULTS.frequency);
        expect(cellCountOf(DIRECTION_DEFAULTS.frequency)).toBe(1002);
        expect(Object.isFrozen(DIRECTION_DEFAULTS)).toBe(true);
        expect(FREQUENCY_OPTIONS[0]).toBe('auto');
        // Every readout from the raw histogram at ν = 10; only the drawn map is smoothed.
        expect(directionReadoutRequest({ ...DIRECTION_DEFAULTS, referenceNumber: 3 })).toEqual({
            referenceNumber: 3, frequency: 10, weight: 'count', frame: 'cartesian', minAmplitudeQuantile: 0, smoothing: 0
        });
        expect(directionMapRequest({ ...DIRECTION_DEFAULTS, referenceNumber: 3 })).toEqual({
            referenceNumber: 3, frequency: 10, weight: 'count', frame: 'cartesian', minAmplitudeQuantile: 0, smoothing: 4
        });
        // Auto is still an absent frequency: the engine resolves it, as for the tool.
        expect(directionReadoutRequest({ ...DIRECTION_DEFAULTS, referenceNumber: 3, frequency: 'auto' }).frequency).toBeNull();
    });

    it('arm A keeps the settings frozen at protocol 2026-10-24.6, not the page\'s current defaults', () => {
        expect(ARM_A_RULES.DIRECTION_SETTINGS).toEqual({ frequency: 'auto', weight: 'count', frame: 'cartesian', minQuantile: 0, smoothing: 0 });
        expect(ARM_A_RULES.DIRECTION_SETTINGS).not.toBe(DIRECTION_DEFAULTS);
        expect(Object.isFrozen(ARM_A_RULES.DIRECTION_SETTINGS)).toBe(true);
        expect(ARM_A_RULES).not.toHaveProperty('FREQUENCY');
        expect(ARM_A_RULES).not.toHaveProperty('SMOOTHING');
        // Its request: the engine's resolution, unsmoothed (the tool's request).
        expect(directionReadoutRequest({ ...ARM_A_RULES.DIRECTION_SETTINGS, referenceNumber: 3 })).toEqual({
            referenceNumber: 3, frequency: null, weight: 'count', frame: 'cartesian', minAmplitudeQuantile: 0, smoothing: 0
        });
    });

    it('analyze_displacement_directions keeps its own defaults: the engine\'s resolution, unsmoothed', async () => {
        const requests = [];
        const request = async (kind, params) => {
            requests.push({ kind, params });
            throw new Error('recorded');
        };
        const tool = buildCopilotTools({ request, structure: null, elements: ['O'] })
            .find((item) => item.name === 'analyze_displacement_directions');
        await expect(tool.run({ site_ref: 3 })).rejects.toThrow('recorded');
        // No frequency (Auto) and no smoothing (the engine's 0): not the page's ν = 10, 4×.
        expect(requests).toEqual([{ kind: 'orientation', params: { geometry: false, referenceNumber: 3 } }]);
        // Its Show in opens the page at Auto, over the page's ν = 10.
        expect(tool.pageParams({ site_ref: 3 })).toEqual({ referenceNumber: 3, frequency: 'auto' });
    });

    it('the readout request never carries the smoothing; the map request differs by the smoothing alone', () => {
        for (const smoothing of [0, 1, 2, 12]) {
            for (const frequency of ['auto', 2, 10, 24]) {
                const settings = { ...DIRECTION_DEFAULTS, referenceNumber: 1, frequency, smoothing };
                expect(directionReadoutRequest(settings)).toEqual(directionReadoutRequest({ ...settings, smoothing: 0 }));
                expect(directionMapRequest(settings)).toEqual({ ...directionReadoutRequest(settings), smoothing });
            }
        }
        expect(directionReadoutRequest({ ...DIRECTION_DEFAULTS, referenceNumber: 1, frequency: 10 }).frequency).toBe(10);
    });

    it('labels Auto with the engine\'s ν for the site on screen, plain Auto before it is known', () => {
        expect(autoResolutionLabel(null)).toBe('Auto');
        expect(autoResolutionLabel({})).toBe('Auto');
        expect(autoResolutionLabel({ recommendedFrequency: 2, frequency: 2 })).toBe('Auto · ν=2 (42)');
        // With a resolution picked, Auto still names its own choice.
        expect(autoResolutionLabel({ recommendedFrequency: 1, frequency: 10 })).toBe('Auto · ν=1 (12)');
        expect(FREQUENCY_OPTIONS.filter((option) => option !== 'auto').map(cellCountOf))
            .toEqual([42, 92, 162, 252, 362, 642, 1002, 1442, 2562, 4002, 5762]);
    });
});
