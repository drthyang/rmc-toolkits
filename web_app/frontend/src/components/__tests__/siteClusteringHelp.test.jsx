// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// The Cluster control (coordinates-only .rmc6f files) is the same control on
// the PCA Ellipsoid and Displacement Directions pages, so its "About site
// clustering" help must say the same thing on both — including what the
// count/copies figure beside a site means and which way to move the distance.

import React from 'react';
import { describe, expect, it, vi } from 'vitest';
import { renderToStaticMarkup } from 'react-dom/server';

vi.mock('../../useSiteCloud', () => ({
    default: () => ({
        sites: {
            reconstructed: true,
            elements: ['Ga'],
            sites: [{ referenceNumber: 1, element: 'Ga', count: 27, copiesPerCell: 27 }],
        },
        sitesError: null,
        loadingSites: false,
        selectedRef: 1,
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

// The popover text of the "About site clustering" InfoBadge, tags stripped.
const clusteringHelp = (Page) => {
    const html = renderToStaticMarkup(React.createElement(Page, { directory: '.', localRun: null }));
    const match = html.match(/aria-label="About site clustering"[^>]*>\?<\/button><span[^>]*role="tooltip"[^>]*>(.*?)<\/span><\/span>/s);
    expect(match).not.toBeNull();
    return match[1].replace(/<[^>]+>/g, ' ').replace(/\s+/g, ' ').trim();
};

describe('About site clustering help', () => {
    const pca = clusteringHelp(PcaKdePage);
    const directions = clusteringHelp(OrientationPage);

    it('explains the count/copies figure and the direction to move the distance', () => {
        expect(pca).toMatch(/the count beside a site \(e\.g\. 27\/27\) is its members against that expected number/);
        expect(pca).toMatch(/Raise the distance to merge over-split sites/);
    });

    it('reads the same on the PCA Ellipsoid and Displacement Directions pages', () => {
        expect(directions).toBe(pca);
    });
});
