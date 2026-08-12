// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

import { describe, it, expect } from 'vitest';

import { moveRatios } from '../moveStats.js';

describe('moveRatios', () => {
    it('divides the move counters by the atom count', () => {
        // The bundled demo model: 1 560 869 generated, 183 014 accepted, 52 000 atoms.
        const ratios = moveRatios({ generated: 1560869, accepted: 183014 }, 52000);
        expect(ratios.generatedPerAtom).toBeCloseTo(30.017, 3);
        expect(ratios.acceptedPerAtom).toBeCloseTo(3.5195, 4);
        expect(ratios.acceptedPerGenerated).toBeCloseTo(0.11725, 5);
    });

    it('keeps the raw counts alongside the ratios', () => {
        expect(moveRatios({ generated: 100, accepted: 25 }, 10)).toMatchObject({
            generated: 100, accepted: 25, generatedPerAtom: 10, acceptedPerAtom: 2.5,
            acceptedPerGenerated: 0.25,
        });
    });

    it('reports a run that accepted nothing as zero, not as missing', () => {
        const ratios = moveRatios({ generated: 500, accepted: 0 }, 100);
        expect(ratios.acceptedPerAtom).toBe(0);
        expect(ratios.acceptedPerGenerated).toBe(0);
    });

    it('leaves the acceptance ratio out when nothing was generated', () => {
        // 0 accepted of 0 generated is undefined, not 0.
        const ratios = moveRatios({ generated: 0, accepted: 0 }, 100);
        expect(ratios.acceptedPerGenerated).toBeUndefined();
        expect(ratios.generatedPerAtom).toBe(0);
    });

    it('handles a partial header, reporting only what it can', () => {
        const ratios = moveRatios({ generated: 400 }, 100);
        expect(ratios.generatedPerAtom).toBe(4);
        expect(ratios.acceptedPerAtom).toBeUndefined();
        expect(ratios.acceptedPerGenerated).toBeUndefined();
    });

    it('returns null when the counters or the atom count are unusable', () => {
        expect(moveRatios(null, 52000)).toBeNull();
        expect(moveRatios(undefined, 52000)).toBeNull();
        expect(moveRatios({ generated: 100 }, 0)).toBeNull();
        expect(moveRatios({ generated: 100 }, undefined)).toBeNull();
        expect(moveRatios({}, 100)).toBeNull();
        expect(moveRatios({ generated: null, accepted: null }, 100)).toBeNull();
    });
});
