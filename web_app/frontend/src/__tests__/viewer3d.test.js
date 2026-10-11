// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// The crystal-direction math behind the Copilot's 3D view tools: how a
// model's spelling of a direction is read, [u v w] ↔ Cartesian in cubic and
// oblique (hexagonal) cells, the low-index name of a view, the screen-up
// convention (down a or b: c up; down c: b up), and the viewer registry.

import { describe, expect, it, vi } from 'vitest';
import {
    angleDeg, cartesianToUvw, describeDirection, directionLabel, directionVector, getViewer, nearestLowIndex,
    normalize, parseDirection, posePlan, registerViewer, registeredViewers, uvwToCartesian, viewUp, waitForViewer
} from '../viewer3d';

const CUBIC = [[10.4, 0, 0], [0, 10.4, 0], [0, 0, 10.4]];
// a = b = 4, c = 6.5, γ = 120°.
const HEXAGONAL = [[4, 0, 0], [-2, 2 * Math.sqrt(3), 0], [0, 0, 6.5]];

const close = (actual, expected, digits = 9) => {
    expect(actual).toHaveLength(expected.length);
    actual.forEach((value, i) => expect(value).toBeCloseTo(expected[i], digits));
};

describe('parseDirection', () => {
    it('reads [uvw] in the spellings models use', () => {
        expect(parseDirection('[111]')).toEqual({ uvw: [1, 1, 1] });
        expect(parseDirection('[1 -1 0]')).toEqual({ uvw: [1, -1, 0] });
        expect(parseDirection('[1-10]')).toEqual({ uvw: [1, -1, 0] });
        expect(parseDirection('1,1,0')).toEqual({ uvw: [1, 1, 0] });
        expect(parseDirection('<110>')).toEqual({ uvw: [1, 1, 0] });
        expect(parseDirection('[1̄10]')).toEqual({ uvw: [-1, 1, 0] });
        expect(parseDirection('[−1 1 0]')).toEqual({ uvw: [-1, 1, 0] });
        expect(parseDirection('[0.71 0.71 1]')).toEqual({ uvw: [0.71, 0.71, 1] });
        expect(parseDirection('[12 3 4]')).toEqual({ uvw: [12, 3, 4] });
        expect(parseDirection([0, 0, 1])).toEqual({ uvw: [0, 0, 1] });
    });

    it('reads axis names', () => {
        expect(parseDirection('c')).toEqual({ axis: 2 });
        expect(parseDirection('a axis')).toEqual({ axis: 0 });
        expect(parseDirection('along b')).toEqual({ axis: 1 });
        expect(parseDirection('c-axis')).toEqual({ axis: 2 });
        expect(parseDirection('PC1')).toEqual({ pc: 0 });
        expect(parseDirection('pc 3')).toEqual({ pc: 2 });
        expect(parseDirection('z')).toEqual({ cartesian: [0, 0, 1] });
    });

    it('refuses what it cannot read', () => {
        expect(parseDirection('[000]')).toBeNull();
        expect(parseDirection('[11]')).toBeNull();
        expect(parseDirection('diagonal')).toBeNull();
        expect(parseDirection('[1 1 x]')).toBeNull();
        expect(parseDirection(null)).toBeNull();
        expect(parseDirection('')).toBeNull();
    });
});

describe('directions in a cell', () => {
    it('maps [u v w] to u·a + v·b + w·c and back, in an oblique cell too', () => {
        close(uvwToCartesian([1, 1, 0], HEXAGONAL), [2, 2 * Math.sqrt(3), 0]);
        close(cartesianToUvw(uvwToCartesian([2, -1, 3], HEXAGONAL), HEXAGONAL), [2, -1, 3]);
        // In a hexagonal cell [110] is not along x + y, and is 60° from a.
        expect(angleDeg(uvwToCartesian([1, 1, 0], HEXAGONAL), HEXAGONAL[0])).toBeCloseTo(60, 9);
        expect(uvwToCartesian([1, 1], CUBIC)).toBeNull();
    });

    it('names a view by its simplest integer direction within 2°, else by its components', () => {
        expect(nearestLowIndex([1, 1, 1], CUBIC)).toEqual({ uvw: [1, 1, 1], offDeg: 0, exact: true });
        expect(nearestLowIndex([2, 2, 2.02], CUBIC)).toMatchObject({ uvw: [1, 1, 1], exact: true });
        expect(nearestLowIndex([0, 0, -3], CUBIC).uvw).toEqual([0, 0, -1]);
        expect(nearestLowIndex(uvwToCartesian([1, 1, 0], HEXAGONAL), HEXAGONAL).uvw).toEqual([1, 1, 0]);
        // 0.43:1:0 has no direction with indices ≤ 4 within 2° ([1 2 0] is 3.3° off).
        expect(nearestLowIndex([0.43, 1, 0], CUBIC)).toEqual({ uvw: [0.43, 1, 0], offDeg: 0, exact: false });
        expect(directionLabel([1, -1, 0])).toBe('[1 -1 0]');
        expect(directionLabel([0.705, -0, 1])).toBe('[0.7 0 1]');
    });

    it('puts c up when looking down a or b, b up looking down c, and c projected otherwise', () => {
        close(viewUp(CUBIC[0], CUBIC), [0, 0, 1]);
        close(viewUp(CUBIC[1], CUBIC), [0, 0, 1]);
        close(viewUp(CUBIC[2], CUBIC), [0, 1, 0]);
        const along111 = normalize([1, 1, 1]);
        const up = viewUp(along111, CUBIC);
        expect(angleDeg(up, along111)).toBeCloseTo(90, 9);
        close(up, normalize([-1, -1, 2]));
        // Looking down hexagonal c: b is up (a and b are both ⟂ c).
        close(viewUp(HEXAGONAL[2], HEXAGONAL), normalize(HEXAGONAL[1]));
        // No cell: z, then y.
        close(viewUp([1, 0, 0]), [0, 0, 1]);
        close(viewUp([0, 0, 1]), [0, 1, 0]);
    });

    it('resolves a parsed direction to the line of sight', () => {
        close(directionVector({ uvw: [1, 1, 1] }, CUBIC), normalize([1, 1, 1]));
        close(directionVector({ axis: 1 }, HEXAGONAL), normalize(HEXAGONAL[1]));
        close(directionVector({ pc: 0 }, CUBIC, [[0, 1, 0], [1, 0, 0], [0, 0, 1]]), [0, 1, 0]);
        expect(directionVector({ pc: 0 }, CUBIC, null)).toBeNull();
        expect(directionVector({ uvw: [1, 0, 0] }, null)).toBeNull();
        close(directionVector({ cartesian: [0, 0, 2] }, null), [0, 0, 1]);
    });

    it('plans a pose on the +direction side of the target', () => {
        const plan = posePlan({ target: [1, 2, 3], direction: [0, 0, 2], distance: 5, up: [0, 1, 0] });
        expect(plan).toEqual({ target: [1, 2, 3], position: [1, 2, 8], direction: [0, 0, 1], up: [0, 1, 0] });
        expect(posePlan({ target: [0, 0, 0], direction: [0, 0, 0], distance: 1 })).toBeNull();
        expect(posePlan({ target: [0, 0, 0], direction: [1, 0, 0], distance: 0 })).toBeNull();
    });

    it('describes a view for the model', () => {
        expect(describeDirection(normalize([1, 1, 1]), normalize([-1, -1, 2]), CUBIC))
            .toEqual({ view_along: '[1 1 1]', screen_up: '[-1 -1 2]' });
        expect(describeDirection([0.43, 1, 0], [0, 0, 1], CUBIC)).toEqual({ view_along: '≈[0.43 1 0]', screen_up: '[0 0 1]' });
        expect(describeDirection([0, 0, 2], null, null)).toEqual({ view_along_cartesian: [0, 0, 1] });
    });
});

describe('the viewer registry', () => {
    it('registers, finds, waits for readiness, and unregisters only its own controller', async () => {
        let ready = false;
        const controller = { ready: () => ready };
        const unregister = registerViewer('test-view', controller);
        expect(getViewer('test-view')).toBe(controller);
        expect(registeredViewers()).toContain('test-view');

        setTimeout(() => { ready = true; }, 30);
        await expect(waitForViewer('test-view', { timeoutMs: 500, pollMs: 10 })).resolves.toBe(controller);
        await expect(waitForViewer('absent', { timeoutMs: 40, pollMs: 10 })).resolves.toBeNull();
        const aborted = new AbortController();
        aborted.abort();
        await expect(waitForViewer('absent', { signal: aborted.signal })).rejects.toMatchObject({ name: 'AbortError' });

        // A remount registers a new controller first; the old cleanup leaves it.
        const next = { ready: vi.fn(() => true) };
        const unregisterNext = registerViewer('test-view', next);
        unregister();
        expect(getViewer('test-view')).toBe(next);
        unregisterNext();
        expect(getViewer('test-view')).toBeNull();
    });
});
