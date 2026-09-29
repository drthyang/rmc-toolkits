// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// planStatRows: how a stat rail lays out stats that may not fit on one line.

import { describe, expect, it } from 'vitest';
import { planStatRows } from '../statRows';

describe('planStatRows', () => {
    it('keeps the one-line row when every stat fits', () => {
        expect(planStatRows([400, 300, 300], 1000)).toEqual({ mode: null, breakA: false, breakB: false });
        expect(planStatRows([900], 900)).toEqual({ mode: null, breakA: false, breakB: false });
    });

    it('uses the aligned grid for a rail without bands that does not fit', () => {
        expect(planStatRows([1200], 900)).toEqual({ mode: 'wrapped', breakA: false, breakB: false });
    });

    it('packs bands greedily, breaking only before a band that does not fit', () => {
        // Model information on a 13" MacBook: cell + counts share a row, the moves drop.
        expect(planStatRows([480, 400, 480], 1250)).toEqual({ mode: 'banded', breakA: false, breakB: true });
        // iPad Pro portrait: the counts drop, the moves fit beside them.
        expect(planStatRows([480, 400, 380], 810)).toEqual({ mode: 'banded', breakA: true, breakB: false });
        // iPhone landscape: one row per band.
        expect(planStatRows([480, 400, 480], 656)).toEqual({ mode: 'banded', breakA: true, breakB: true });
    });

    it('handles a single band start', () => {
        expect(planStatRows([700, 500], 1000)).toEqual({ mode: 'banded', breakA: true, breakB: false });
    });
});
