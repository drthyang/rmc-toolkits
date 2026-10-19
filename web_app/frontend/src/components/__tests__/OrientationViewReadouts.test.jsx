// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// Plan item 2.5: the readouts under the direction sphere are the raw
// histogram's (`readout`, smoothing 0, what analyze_displacement_directions
// reads), never the drawn map's. The engines pick the peak cell on the
// smoothed map, so with smoothing on the drawn map's peak can sit elsewhere,
// at another σ; the colorbar and the legend still describe the drawn map.

import { describe, expect, it } from 'vitest';
import { renderToStaticMarkup } from 'react-dom/server';
import OrientationView from '../OrientationView';

const common = {
    cellCount: 42,
    orientationAnisotropy: 0.41,
    orientationAnisotropyNull: 0.07,
    orientationAnisotropySignificance: 9.7,
    antipodalAsymmetry: 0.12,
    antipodalAsymmetryNull: 0.1,
    antipodalAsymmetryNullSd: 0.02,
    antipodalAsymmetryZ: 1,
    antipodalAsymmetrySignificant: false,
    mapSignificance: 7.2,
    mapChiSquare: 210,
    mapDegreesOfFreedom: 41,
    mapNullSd: 9,
    weight: 'count',
    browserOrientation: true
};
// A raw histogram and the map drawn from it with two smoothing passes
// (synthetic numbers; ν = 2, whose cells reach 19.3–20.9° from their centres).
const readout = {
    ...common, frequency: 2, smoothing: 0, vmax: 2.5, peakEnhancement: 2.5, peakDirection: [0.6, -0.8, 0], peakCellRadiusDeg: 20.905,
    peakTieCount: 1, peakSignificance: 3.4
};
const map = {
    ...common, frequency: 2, smoothing: 2, vmax: 1.5, peakEnhancement: 1.5, peakDirection: [-0.8, 0, 0.6], peakCellRadiusDeg: 19.348,
    peakTieCount: 3, peakSignificance: -2.1
};

const render = (props) => renderToStaticMarkup(
    <OrientationView
        loadedRef={3}
        selectedRef={3}
        selectedEllipsoid={null}
        unitCell={null}
        frame="cartesian"
        onFrameChange={() => {}}
        colormap="viridis"
        contrast={1}
        relief={0.5}
        showOutline
        showAxes
        {...props}
    />
);
const strip = (html) => html.slice(html.indexOf('ui-inline-stats'), html.indexOf('ui-legend'));

describe('the direction sphere\'s readouts are the raw histogram\'s', () => {
    it('the peak, its direction and σ come from the readout, the colorbar from the drawn map', () => {
        const html = render({ result: map, readout });
        expect(strip(html)).toContain('2.50×');
        expect(strip(html)).toContain('[0.60, -0.80, 0.00]');
        expect(strip(html)).toContain('3.4σ');
        // The direction with its resolution: the peak cell's radius (plan item 2.5).
        expect(strip(html)).toContain(' ± 20.9°, the cell radius at ν = 2');
        // Nothing of the smoothed map's own peak.
        expect(strip(html)).not.toContain('1.50×');
        expect(strip(html)).not.toContain('19.3°');
        expect(strip(html)).not.toContain('equal cells');
        expect(strip(html)).not.toContain('not significant');
        // The colorbar and the legend describe what is drawn.
        expect(html).toContain('1.5× isotropic');
        expect(html).toContain('smoothed 2×');
    });

    it('unsmoothed, the readout is the map and reads the same', () => {
        expect(strip(render({ result: readout, readout }))).toBe(strip(render({ result: map, readout })));
    });

    it('no readout, no readouts: never the drawn map\'s numbers in their place', () => {
        const html = render({ result: map, readout: null });
        expect(html).toContain('orient-colorbar-row');
        expect(html).not.toContain('ui-inline-stats');
    });
});
