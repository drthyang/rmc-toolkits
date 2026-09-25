// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// One label rule for a PCA/orientation site on both pages that show sites.

import { describe, expect, it } from 'vitest';
import { siteLabel } from '../siteLabel';

describe('siteLabel', () => {
    it('names a pure site by its element', () => {
        expect(siteLabel({ element: 'Se', mixed: false, elementCounts: { Se: 64 } })).toBe('Se');
        expect(siteLabel({ element: 'Nb' })).toBe('Nb');
        expect(siteLabel(null)).toBe('');
    });

    it('names a mixed site by its composition, majority first and ties by name', () => {
        expect(siteLabel({ element: 'Ga', mixed: true, elementCounts: { In: 16, Ga: 48 } })).toBe('Ga0.75In0.25');
        expect(siteLabel({ element: 'Co', mixed: true, elementCounts: { Fe: 32, Co: 32 } })).toBe('Co0.50Fe0.50');
    });
});
