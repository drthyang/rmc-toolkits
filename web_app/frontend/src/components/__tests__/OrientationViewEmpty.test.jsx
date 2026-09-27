// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// With no site (no run yet, or no site picked) the Displacement Directions
// sphere card must not be a blank canvas: it carries the same one-line empty
// state as the PCA Ellipsoid table card, "No site selected.".

import { describe, expect, it } from 'vitest';
import { renderToStaticMarkup } from 'react-dom/server';
import OrientationView from '../OrientationView';

const render = () => renderToStaticMarkup(
    <OrientationView
        requestPca={async () => null}
        ready={false}
        selectedRef={null}
        selectedEllipsoid={null}
        clusterThreshold={1}
        unitCell={null}
        frequency="auto"
        weight="count"
        frame="cartesian"
        onFrameChange={() => {}}
        smoothing={0}
        minQuantile={0}
        colormap="viridis"
        contrast={1}
        relief={0.5}
        showOutline
        showAxes
    />
);

describe('OrientationView empty state', () => {
    it('says "No site selected." in the sphere card, as the PCA page does', () => {
        const html = render();
        const sphereCard = html.slice(html.indexOf('orient-main-panel'));
        expect(sphereCard).toContain('<p class="ui-card__caption">No site selected.</p>');
        // Once, and not in the axis-views card.
        expect(html.split('No site selected.')).toHaveLength(2);
    });

    it('shows no readouts, loading badge or error without a result', () => {
        const html = render();
        expect(html).not.toContain('Computing…');
        expect(html).not.toContain('ui-overlay-badge');
        expect(html).not.toContain('orient-colorbar-row');
    });
});
