// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// While the site table loads, the shared Site-ellipsoids panel carries its
// "Loading…" chip on BOTH pages that mount it — the PCA Ellipsoid page used to
// leave the prop out, so its picker sat empty with no cue.

import React from 'react';
import { describe, expect, it, vi } from 'vitest';
import { renderToStaticMarkup } from 'react-dom/server';

const state = vi.hoisted(() => ({ loadingSites: true }));

vi.mock('../../useSiteCloud', () => ({
    default: () => ({
        sites: null,
        sitesError: null,
        loadingSites: state.loadingSites,
        selectedRef: null,
        setSelectedRef: () => {},
        selectedEllipsoid: null,
        requestPca: async () => null,
        localFile: null,
        rmc6fText: null,
        ready: false,
        unitCell: null,
        datasetKey: null,
    }),
}));

const { default: PcaKdePage } = await import('../PcaKdePage');
const { default: OrientationPage } = await import('../OrientationPage');

// The markup from the site panel's root to the end of its header bar.
const sitePanelHeader = (html) => {
    const start = html.indexOf('pca-unitcell-panel');
    expect(start).toBeGreaterThan(-1);
    const panel = html.slice(start);
    return panel.slice(0, panel.indexOf('pca-structure'));
};

const LOADING_CHIP = '<span class="ui-card__meta">Loading…</span>';

const PAGES = { 'PCA Ellipsoid': PcaKdePage, 'Displacement Directions': OrientationPage };
const render = (name) => renderToStaticMarkup(React.createElement(PAGES[name], { directory: '.', localRun: null }));

describe.each(Object.keys(PAGES))('%s site panel', (name) => {
    it('shows the Loading… chip while the sites load', () => {
        state.loadingSites = true;
        expect(sitePanelHeader(render(name))).toContain(LOADING_CHIP);
    });

    it('drops it once they have loaded', () => {
        state.loadingSites = false;
        expect(sitePanelHeader(render(name))).not.toContain(LOADING_CHIP);
    });
});
