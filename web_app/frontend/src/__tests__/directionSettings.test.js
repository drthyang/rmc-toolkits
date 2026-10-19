// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// Plan item 2.5: the Displacement Directions page's settings, one module for
// the page and the validation's arm A. The page opens at the engine's
// resolution with the tool's analysis settings, and its readouts never carry
// the smoothing (a display setting).

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

describe('the Displacement Directions page settings (plan item 2.5)', () => {
    it('open at Auto, with the analysis settings analyze_displacement_directions uses', () => {
        expect(DIRECTION_DEFAULTS).toEqual({ frequency: 'auto', weight: 'count', frame: 'cartesian', minQuantile: 0, smoothing: 2 });
        expect(Object.isFrozen(DIRECTION_DEFAULTS)).toBe(true);
        expect(FREQUENCY_OPTIONS[0]).toBe('auto');
        // No frequency: the engine resolves it to its recommended one, as for the tool.
        expect(directionReadoutRequest({ ...DIRECTION_DEFAULTS, referenceNumber: 3 })).toEqual({
            referenceNumber: 3, frequency: null, weight: 'count', frame: 'cartesian', minAmplitudeQuantile: 0, smoothing: 0
        });
    });

    it('arm A reads the page\'s defaults themselves (one definition, never a copied number)', () => {
        expect(ARM_A_RULES.DIRECTION_SETTINGS).toBe(DIRECTION_DEFAULTS);
        expect(ARM_A_RULES).not.toHaveProperty('FREQUENCY');
        expect(ARM_A_RULES).not.toHaveProperty('SMOOTHING');
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
