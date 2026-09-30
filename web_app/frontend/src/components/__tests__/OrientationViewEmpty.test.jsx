// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// With no site (no run yet, or no site picked) the Displacement Directions
// sphere card must not be a blank canvas: its stage carries the PCA Ellipsoid
// page's wording, "Loading sites…" while the sites load, then "No site selected.".

import { describe, expect, it } from 'vitest';
import { renderToStaticMarkup } from 'react-dom/server';
import OrientationView from '../OrientationView';

const render = (extra = {}) => renderToStaticMarkup(
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
        {...extra}
    />
);

describe('OrientationView empty state', () => {
    it('says "No site selected." in the sphere stage, as the PCA page does', () => {
        const html = render();
        const stage = html.slice(html.indexOf('orient-canvas'));
        expect(stage).toContain('<div class="ui-overlay-badge">No site selected.</div>');
        // Once, and not in the axis-views card.
        expect(html.split('No site selected.')).toHaveLength(2);
    });

    it('says "Loading sites…" while the sites load', () => {
        const html = render({ loadingSites: true });
        expect(html).toContain('<div class="ui-overlay-badge">Loading sites…</div>');
        expect(html).not.toContain('No site selected.');
    });

    it('shows no readouts, computing badge or error without a result', () => {
        const html = render();
        expect(html).not.toContain('Computing…');
        expect(html).not.toContain('is-error');
        expect(html).not.toContain('orient-colorbar-row');
    });
});
