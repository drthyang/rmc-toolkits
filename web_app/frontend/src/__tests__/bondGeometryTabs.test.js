// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

import { describe, expect, it } from 'vitest';
import { grTabForFocus } from '../bondGeometryTabs';

describe('grTabForFocus: a Copilot highlight brings its tab forward', () => {
    it('bond-windows asks for the g(r) tab, bond-angles for the angle view', () => {
        expect(grTabForFocus({ panel: 'primary' })).toBe(true);
        expect(grTabForFocus({ panel: 'secondary' })).toBe(true);
        expect(grTabForFocus({ card: 'primary' })).toBe(false);
        expect(grTabForFocus({ card: 'secondary' })).toBe(false);
    });

    it('both named: the primary one wins; two secondaries change nothing', () => {
        expect(grTabForFocus({ card: 'secondary', panel: 'primary' })).toBe(true);
        expect(grTabForFocus({ card: 'primary', panel: 'secondary' })).toBe(false);
        expect(grTabForFocus({ card: 'secondary', panel: 'secondary' })).toBeNull();
    });

    it('a running analysis, a hovered chip or no focus changes nothing', () => {
        expect(grTabForFocus({ card: 'running' })).toBeNull();
        expect(grTabForFocus({ panel: 'preview' })).toBeNull();
        expect(grTabForFocus({})).toBeNull();
        expect(grTabForFocus()).toBeNull();
    });
});
