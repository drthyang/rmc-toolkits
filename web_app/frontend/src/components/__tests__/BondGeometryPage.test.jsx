// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// The Bond Geometry page's in-app help must describe what the page draws.

import { describe, expect, it } from 'vitest';
import { renderToStaticMarkup } from 'react-dom/server';
import BondGeometryPage from '../BondGeometryPage';

const html = renderToStaticMarkup(<BondGeometryPage directory="." localRun={null} />);
const text = html.replace(/<[^>]+>/g, '').replace(/&#x27;/g, "'").replace(/\s+/g, ' ');

describe('Bond Geometry help text (triplets.physics.10 / physics.23 / numerics.33)', () => {
    it('the window helper marks the window with dashed guides; nothing is shaded', () => {
        expect(text).not.toMatch(/shades/);
        expect(text).toMatch(/dashed guides/);
    });

    it('the second partial depends on the bond types, not on the Distinct B–C switch', () => {
        expect(text).toMatch(/second partial is drawn whenever A–B and B–C are different pair types/);
        expect(text).not.toMatch(/With distinct B–C on, the B–C partial is plotted/);
    });

    it('does not present sin-corrected as RMCProfile\'s own normalization', () => {
        expect(text).not.toMatch(/This is RMCProfile's sinth view/);
        expect(text).toMatch(/norm\/sin\(theta\)/);
    });
});
