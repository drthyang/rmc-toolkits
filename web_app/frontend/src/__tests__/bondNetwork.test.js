// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// The folded cell's bonds (bondNetwork.js): the linked-cell search finds
// exactly the bonds the panel's earlier all-pairs loop drew, the draw cap
// keeps the bonds nearest the centre and says so, and the selection's colours
// and fading follow the bond sets.

import { describe, expect, it } from 'vitest';
import {
    bondRadius, buildBondNetwork, cloudOpacity, CLOUD_OPACITY, CLOUD_OPACITY_FADED, findSiteBonds, MAX_DRAWN_BONDS,
} from '../bondNetwork';
import { BOND_COLORS } from '../plotPalette';

// The panel's earlier search (FoldedCellPanel before the linked cells): every
// A site against every B site and 27 images.
const bruteForce = (sites, unit, set) => {
    const [elementA, elementB] = set.elements;
    const [rmin, rmax] = set.window;
    const cart = (f) => [0, 1, 2].map((axis) => f[0] * unit[0][axis] + f[1] * unit[1][axis] + f[2] * unit[2][axis]);
    const out = [];
    sites.forEach((from) => {
        if (from.element !== elementA) return;
        sites.forEach((to) => {
            if (to.element !== elementB) return;
            for (let mx = -1; mx <= 1; mx += 1) {
                for (let my = -1; my <= 1; my += 1) {
                    for (let mz = -1; mz <= 1; mz += 1) {
                        const inCell = mx === 0 && my === 0 && mz === 0;
                        if (inCell && to.referenceNumber === from.referenceNumber) continue;
                        if (elementA === elementB && inCell && to.referenceNumber < from.referenceNumber) continue;
                        const delta = cart([
                            to.siteFractional[0] + mx - from.siteFractional[0],
                            to.siteFractional[1] + my - from.siteFractional[1],
                            to.siteFractional[2] + mz - from.siteFractional[2],
                        ]);
                        const length = Math.hypot(...delta);
                        if (length < rmin || length > rmax || length === 0) continue;
                        const start = cart(from.siteFractional);
                        out.push([...start, ...start.map((value, axis) => value + delta[axis])]);
                    }
                }
            }
        });
    });
    return out;
};

const rows = (segments) => {
    const out = [];
    for (let i = 0; i < segments.length; i += 6) out.push([...segments.subarray(i, i + 6)]);
    return out;
};
const key = (row) => row.map((value) => value.toFixed(6)).join(',');
const sameBonds = (a, b) => expect(a.map(key).sort()).toEqual(b.map(key).sort());

// A deterministic pseudo-random generator (mulberry32).
const random = (seed) => () => {
    let t = (seed += 0x6d2b79f5);
    t = Math.imul(t ^ (t >>> 15), t | 1);
    t ^= t + Math.imul(t ^ (t >>> 7), t | 61);
    return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
};

// A rock-salt-like cell: n×n×n sub-cells of spacing a, A and B alternating.
const rockSalt = (n, a = 2.5) => {
    const sites = [];
    let ref = 1;
    for (let i = 0; i < n; i += 1) {
        for (let j = 0; j < n; j += 1) {
            for (let k = 0; k < n; k += 1) {
                sites.push({ referenceNumber: ref++, element: (i + j + k) % 2 ? 'B' : 'A', siteFractional: [i / n, j / n, k / n] });
            }
        }
    }
    const L = n * a;
    return { sites, unit: [[L, 0, 0], [0, L, 0], [0, 0, L]] };
};

describe('findSiteBonds: the linked-cell search draws the same bonds as the all-pairs loop', () => {
    it('on a triclinic cell, for different and same elements, windows short and past the cell', () => {
        const next = random(7);
        const unit = [[6.1, 0, 0], [1.3, 5.4, 0], [-0.8, 1.1, 7.2]];
        const sites = Array.from({ length: 40 }, (_, index) => ({
            referenceNumber: index + 1,
            element: ['Ga', 'Ta', 'Se'][index % 3],
            siteFractional: [next(), next(), next()],
        }));
        for (const set of [
            { elements: ['Ta', 'Se'], window: [2.0, 3.0] },
            { elements: ['Se', 'Se'], window: [1.5, 3.4] },
            { elements: ['Ga', 'Ta'], window: [0, 8.5] },
        ]) {
            const fast = rows(findSiteBonds(sites, unit, set));
            const slow = bruteForce(sites, unit, set);
            expect(slow.length).toBeGreaterThan(0);
            sameBonds(fast, slow);
        }
    });

    it('an element missing or an empty window draws nothing', () => {
        const { sites, unit } = rockSalt(2);
        expect(findSiteBonds(sites, unit, { elements: ['A', 'C'], window: [2, 3] })).toHaveLength(0);
        expect(findSiteBonds(sites, unit, { elements: ['A', 'B'], window: [3, 2] })).toHaveLength(0);
        expect(findSiteBonds(sites, unit, { elements: ['A', 'B'], window: [1, 2] })).toHaveLength(0);
    });

    it('a rock-salt cell: six bonds per A site', () => {
        const { sites, unit } = rockSalt(4);
        const bonds = findSiteBonds(sites, unit, { elements: ['A', 'B'], window: [2.4, 2.6] });
        expect(bonds.length / 6).toBe(32 * 6);
    });
});

describe('buildBondNetwork: what the 3D view draws', () => {
    it('keeps each set\'s pair colour and draws every bond under the cap', () => {
        const { sites, unit } = rockSalt(4);
        const network = buildBondNetwork(sites, unit, [
            { elements: ['A', 'B'], window: [2.4, 2.6], color: BOND_COLORS.ab },
            { elements: ['B', 'B'], window: [3.4, 3.6], color: BOND_COLORS.bc },
        ]);
        expect(network.limited).toBe(false);
        expect(network.sets.map((set) => set.color)).toEqual([BOND_COLORS.ab, BOND_COLORS.bc]);
        expect(network.shown).toBe(network.total);
        expect(network.total).toBe(32 * 6 + network.sets[1].segments.length / 6);
    });

    it('a large box: the cap keeps the bonds nearest the centre, across sets, and counts them', () => {
        const { sites, unit } = rockSalt(20); // 8000 sites, a 50 Å cell
        const center = [25, 25, 25];
        const sets = [
            { elements: ['A', 'B'], window: [2.4, 2.6], color: BOND_COLORS.ab },
            { elements: ['A', 'A'], window: [3.4, 3.6], color: BOND_COLORS.bc },
        ];
        const network = buildBondNetwork(sites, unit, sets, { center, maxBonds: 500 });
        expect(network.limited).toBe(true);
        const each = sets.map((set) => findSiteBonds(sites, unit, set).length / 6);
        expect(each[0]).toBe(4000 * 6);
        expect(network.total).toBe(each[0] + each[1]);
        expect(network.shown).toBe(500);
        const drawn = network.sets.reduce((sum, set) => sum + set.segments.length / 6, 0);
        expect(drawn).toBe(500);
        // Every kept bond is nearer the centre than every dropped one.
        const radius = (s, b) => Math.hypot(...[0, 1, 2].map((axis) => (s[b * 6 + axis] + s[b * 6 + 3 + axis]) / 2 - center[axis]));
        const kept = network.sets.flatMap((set) => Array.from({ length: set.segments.length / 6 }, (_, b) => radius(set.segments, b)));
        expect(Math.max(...kept)).toBeLessThan(12);
        expect(network.sets[0].segments.length).toBeGreaterThan(0);
        expect(network.sets[1].segments.length).toBeGreaterThan(0);
    });

    it('the default cap is the panel\'s', () => {
        expect(MAX_DRAWN_BONDS).toBe(1500);
        const { sites, unit } = rockSalt(20);
        const network = buildBondNetwork(sites, unit, [{ elements: ['A', 'B'], window: [2.4, 2.6] }], { center: [25, 25, 25] });
        expect([network.total, network.shown, network.limited]).toEqual([24000, 1500, true]);
    });

    it('no sets, no bonds', () => {
        expect(buildBondNetwork([], [[1, 0, 0], [0, 1, 0], [0, 0, 1]], null)).toEqual({ sets: [], shown: 0, total: 0, limited: false, center: [0, 0, 0], radiusSq: null });
    });
});

describe('the selection stands out', () => {
    it('the tubes are thick enough to see: a fraction of the bond, a floor tied to the cell', () => {
        expect(bondRadius([2, 3], 10.4)).toBeCloseTo(0.0875, 6);
        expect(bondRadius([2, 3], 100)).toBeCloseTo(0.25, 6);
    });

    it('the clouds of elements outside the drawn bonds fade', () => {
        const sets = [{ elements: ['Ta', 'Se'] }];
        expect(cloudOpacity('Ta', sets)).toBe(CLOUD_OPACITY);
        expect(cloudOpacity('Se', sets)).toBe(CLOUD_OPACITY);
        expect(cloudOpacity('Ga', sets)).toBe(CLOUD_OPACITY_FADED);
        // Before Compute (no bonds) nothing fades.
        expect(cloudOpacity('Ga', null)).toBe(CLOUD_OPACITY);
        expect(cloudOpacity('Ga', [])).toBe(CLOUD_OPACITY);
    });
});
