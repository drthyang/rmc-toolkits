// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// Plan item 2.5: the direction resolution is a measurement limit. The peak
// direction is the peak cell's centre, so analyze_displacement_directions
// gives it with the cell's circumscribed radius (the engine's
// peakCellRadiusDeg) and, where the cell is wider than LOW_INDEX_CONE_DEG,
// says the direction is not resolved finer than that and names no low-index
// family. Known answers on a planted [0 0 1] split site (synthetic, through
// the browser engine and the tool's summary).

import { describe, expect, it } from 'vitest';
import { orientationHistogram } from '../workers/orientation';
import { directionResolved, LOW_INDEX_CONE_DEG, summarizeOrientation } from '../copilotTools';
import { mulberry32 } from '../../validation/cases';

const CUBIC = [[4, 0, 0], [0, 4, 0], [0, 0, 4]];
const normal = (rng) => {
    let u = 0;
    while (u === 0) u = rng();
    return Math.sqrt(-2 * Math.log(u)) * Math.cos(2 * Math.PI * rng());
};
// A site split ±delta along [0 0 1] over Gaussian thermal motion (Å).
const splitSite = (n, seed, { delta = 0.15, sigma = 0.05 } = {}) => {
    const rng = mulberry32(seed);
    return Array.from({ length: n }, () => {
        const well = rng() < 0.5 ? delta : -delta;
        return [sigma * normal(rng), sigma * normal(rng), well + sigma * normal(rng)];
    });
};
// The angle (degrees) between a direction and the [0 0 1] axis, either sense.
const offAxis = (direction) => (Math.acos(Math.min(1, Math.abs(direction[2]) / Math.hypot(...direction))) * 180) / Math.PI;
const preferred = (checks) => checks.find((check) => check.text.startsWith('Preferred direction'))?.text;

describe('the peak direction with its cell radius (plan item 2.5)', () => {
    it('216 copies of a planted [0 0 1] split site: ν = 1, significant, but not resolved finer than 37.4°', () => {
        const payload = orientationHistogram(splitSite(216, 41), { geometry: false });
        // Auto: fewer than 504 points give the 12-cell dodecahedron, whose
        // every cell reaches arccos(sqrt((5 + 2√5)/15)) = 37.377° from its centre.
        expect(payload.frequency).toBe(1);
        expect(payload.peakCellRadiusDeg).toBeCloseTo(37.3774, 4);

        const { data, checks } = summarizeOrientation(payload, CUBIC);
        expect(data.peak.radius_deg).toBe(37.4);
        expect(data.peak.significance_sigma).toBeGreaterThanOrEqual(3);
        // The lobe is real, but the cell centre it names sits 31.7° off the
        // planted axis (a cube axis is a cell boundary at ν = 1), within its radius.
        expect(offAxis(data.peak.direction_cartesian)).toBeCloseTo(31.7, 0);
        expect(offAxis(data.peak.direction_cartesian)).toBeLessThanOrEqual(data.peak.radius_deg);
        expect(directionResolved(data.peak.radius_deg)).toBe(false);

        const text = preferred(checks);
        expect(text).toContain(`toward [${data.peak.direction_uvw.join(' ')}] ± 37.4°, the cell radius at ν = 1.`);
        expect(text).toContain('The direction is not resolved finer than 37.4° (over 15°), so it names no <100>, <110> or <111> family.');
        // ...and names none: the only families in the text are the ones it declines.
        expect(text.match(/<\d{3}>/g)).toEqual(['<100>', '<110>', '<111>']);
    });

    it('4000 copies of the same site: Auto is ν = 5 and the direction is resolved, ± its cell radius, near [0 0 1]', () => {
        const payload = orientationHistogram(splitSite(4000, 42), { geometry: false });
        expect(payload.frequency).toBe(5);
        expect(payload.peakCellRadiusDeg).toBeLessThan(LOW_INDEX_CONE_DEG);

        const { data, checks } = summarizeOrientation(payload, CUBIC);
        expect(directionResolved(data.peak.radius_deg)).toBe(true);
        expect(offAxis(data.peak.direction_cartesian)).toBeLessThanOrEqual(data.peak.radius_deg);
        const text = preferred(checks);
        expect(text).toBe(`Preferred direction: [${data.peak.direction_uvw.join(' ')}] ± ${data.peak.radius_deg}°, the cell radius at ν = 5; `
            + `the peak cell is ${data.peak.enhancement}× isotropic at ${data.peak.significance_sigma}σ (corrected over all cells).`);
    });

    it('the radius is in the data whatever the significance; no check claims a direction without it', () => {
        const isotropic = Array.from({ length: 600 }, (_, index) => {
            const rng = mulberry32(index + 1);
            return [normal(rng), normal(rng), normal(rng)];
        });
        const { data, checks } = summarizeOrientation(orientationHistogram(isotropic, { geometry: false }), CUBIC);
        expect(data.frequency).toBe(2);
        expect(data.peak.radius_deg).toBeGreaterThan(19);
        expect(data.peak.significance_sigma).toBeLessThan(3);
        expect(preferred(checks)).toBeUndefined();
    });
});
