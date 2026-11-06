// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// A pointer readout sits right of the point, and left of it past the stage's
// middle, so a tap near a phone's right edge keeps it on the card.

import { describe, expect, it } from 'vitest';
import { tooltipPlace } from '../components/tooltipPlace';

describe('tooltipPlace', () => {
    it('sits right of and below the point in the left half', () => {
        expect(tooltipPlace({ x: 100, y: 40, w: 370 })).toEqual({ left: 114, top: 52 });
    });

    it('sits left of the point past the middle, measured from the right edge', () => {
        expect(tooltipPlace({ x: 300, y: 40, w: 370 })).toEqual({ right: 84, top: 52 });
    });

    it('stays right of the point when the stage width is unknown', () => {
        expect(tooltipPlace({ x: 900, y: 5 })).toEqual({ left: 914, top: 17 });
    });
});
