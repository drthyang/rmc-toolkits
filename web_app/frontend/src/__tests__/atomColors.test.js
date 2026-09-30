// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// One element -> color map per model, the same on every page, with similar
// colors separated automatically: every pair clears MIN_COLOR_DISTANCE in
// OKLab for normal vision and MIN_CVD_COLOR_DISTANCE under simulated
// deuteranopia and protanopia; the more abundant species keeps its table
// color; a model whose colors are already distinct is left alone.

import { readFileSync } from 'node:fs';
import { describe, expect, it } from 'vitest';
import {
    ELEMENT_COLORS,
    MIN_COLOR_DISTANCE,
    MIN_CVD_COLOR_DISTANCE,
    buildElementColors,
    colorDistances,
    elementColorNote,
    resolveElementColors,
    speciesCounts
} from '../atomColors';
import { structureFromRmc6f } from '../browserData';
import { siteDisplacementsFromRmc6f, siteEllipsoids, siteSpecies } from '../workers/pcaKde.js';

const worstPair = (colors) => {
    const names = Object.keys(colors);
    let worst = { normal: Infinity, deuteranopia: Infinity, protanopia: Infinity };
    for (let i = 0; i < names.length; i += 1) {
        for (let j = i + 1; j < names.length; j += 1) {
            const d = colorDistances(colors[names[i]], colors[names[j]]);
            worst = {
                normal: Math.min(worst.normal, d.normal),
                deuteranopia: Math.min(worst.deuteranopia, d.deuteranopia),
                protanopia: Math.min(worst.protanopia, d.protanopia)
            };
        }
    }
    return worst;
};

const MODELS = {
    GaTa4Se8: { Ga: 4000, Ta: 16000, Se: 32000 },
    GaNb4Se8: { Ga: 4000, Nb: 16000, Se: 32000 },
    Mn3Sn: { Mn: 3, Sn: 1 },
    FeCoSn: { Fe: 1, Co: 1, Sn: 1 },
    Fe2O3: { Fe: 2, O: 3 },
    BaTiO3: { Ba: 1, Ti: 1, O: 3 },
    TiVO4: { Ti: 1, V: 1, O: 4 },
    CoCrFeMnNi: { Co: 1, Cr: 1, Fe: 1, Mn: 1, Ni: 1 },
    'no table colors': { Lu: 2, Hf: 1, Re: 1, Os: 1, Ir: 1 }
};

describe('resolveElementColors', () => {
    it.each(Object.entries(MODELS))('%s: every pair is far enough apart, in normal and colour-blind vision', (name, counts) => {
        const { colors } = resolveElementColors(Object.keys(counts), counts);
        const worst = worstPair(colors);
        expect(worst.normal).toBeGreaterThanOrEqual(MIN_COLOR_DISTANCE);
        expect(worst.deuteranopia).toBeGreaterThanOrEqual(MIN_CVD_COLOR_DISTANCE);
        expect(worst.protanopia).toBeGreaterThanOrEqual(MIN_CVD_COLOR_DISTANCE);
    });

    it.each(['GaTa4Se8', 'GaNb4Se8', 'BaTiO3'])('%s: already distinct, so every table color is kept', (name) => {
        const counts = MODELS[name];
        const { colors, adjusted } = resolveElementColors(Object.keys(counts), counts);
        expect(adjusted).toEqual([]);
        Object.keys(counts).forEach((element) => expect(colors[element]).toBe(ELEMENT_COLORS[element]));
    });

    it('moves the minority species of a close pair and keeps the majority one', () => {
        // Fe #e06633 and O #e6443b are 0.06 apart. O is the majority in Fe2O3.
        const oxide = resolveElementColors(['Fe', 'O'], { Fe: 2, O: 3 });
        expect(oxide.colors.O).toBe(ELEMENT_COLORS.O);
        expect(oxide.colors.Fe).not.toBe(ELEMENT_COLORS.Fe);
        expect(oxide.adjusted).toEqual([{ element: 'Fe', from: ELEMENT_COLORS.Fe, to: oxide.colors.Fe, near: 'O' }]);
        // An iron-rich model keeps Fe's color and moves O instead.
        const metal = resolveElementColors(['Fe', 'O'], { Fe: 9, O: 1 });
        expect(metal.colors.Fe).toBe(ELEMENT_COLORS.Fe);
        expect(metal.colors.O).not.toBe(ELEMENT_COLORS.O);
    });

    it('moves a color only as far as it needs to', () => {
        // Ti/V greys: V moves to the nearest color that clears Ti, not across the wheel.
        const { colors } = resolveElementColors(['Ti', 'V'], { Ti: 1, V: 1 });
        const shift = colorDistances(colors.V, ELEMENT_COLORS.V).normal;
        expect(shift).toBeGreaterThan(0);
        expect(shift).toBeLessThan(0.2);
    });

    it('ties go to table colors, then names; atom order and duplicates do not matter', () => {
        const a = resolveElementColors(['Fe', 'Co', 'Sn'], { Fe: 1, Co: 1, Sn: 1 });
        const b = resolveElementColors(['Sn', 'Co', 'Fe', 'Co', 'Sn'], { Sn: 1, Fe: 1, Co: 1 });
        expect(b.colors).toEqual(a.colors);
        // Equal shares: Co is placed before Fe by name, so Fe is the one that moves.
        expect(a.colors.Co).toBe(ELEMENT_COLORS.Co);
        expect(a.adjusted.map((move) => move.element)).toEqual(['Fe']);
    });

    it('lists elements in name order, whatever the placement order', () => {
        const { colors } = resolveElementColors(['Se', 'Ta', 'Ga'], MODELS.GaTa4Se8);
        expect(Object.keys(colors)).toEqual(['Ga', 'Se', 'Ta']);
    });

    it('never picks a near-white or near-black replacement', () => {
        const { colors, adjusted } = resolveElementColors(Object.keys(MODELS.CoCrFeMnNi), MODELS.CoCrFeMnNi);
        expect(adjusted.length).toBeGreaterThan(0);
        adjusted.forEach(({ element }) => {
            expect(colorDistances(colors[element], '#ffffff').normal).toBeGreaterThan(0.15);
            expect(colorDistances(colors[element], '#000000').normal).toBeGreaterThan(0.45);
        });
    });
});

describe('buildElementColors', () => {
    it('agrees with resolveElementColors and notes each move for the legend hover', () => {
        const map = buildElementColors(['Fe', 'O'], { Fe: 2, O: 3 });
        expect({ ...map }).toEqual(resolveElementColors(['Fe', 'O'], { Fe: 2, O: 3 }).colors);
        expect(Object.keys(map)).toEqual(['Fe', 'O']);
        expect(elementColorNote(map, 'Fe')).toContain("too close to O's");
        expect(elementColorNote(map, 'O')).toBeUndefined();
    });

    it('gives the same map when two payloads differ by a few skipped atoms', () => {
        const full = buildElementColors(['Co', 'Fe', 'Sn'], { Co: 1000, Fe: 1000, Sn: 1000 });
        const skipped = buildElementColors(['Co', 'Fe', 'Sn'], { Co: 998, Fe: 1000, Sn: 1000 });
        expect(skipped).toEqual(full);
    });
});

describe('speciesCounts', () => {
    it('sums every site composition, minority species of mixed sites included', () => {
        expect(speciesCounts([
            { element: 'Nb', count: 10, elementCounts: { Nb: 7, Ta: 3 }, mixed: true },
            { element: 'Nb', count: 10, elementCounts: { Nb: 10 } },
            { element: 'Se', count: 40 }
        ])).toEqual({ Nb: 17, Ta: 3, Se: 40 });
        expect(speciesCounts(null)).toEqual({});
    });
});

describe('the demo run: one map on every page', () => {
    const text = readFileSync(new URL('../../public/demo/GTS_250K.rmc6f', import.meta.url), 'utf8');

    it('Atomic Density (structure) and the site pages (PCA, Directions, Bond Geometry) agree', () => {
        const structure = structureFromRmc6f({ path: 'Demo/GTS_250K.rmc6f', text });
        const parsed = siteDisplacementsFromRmc6f(text);
        const sites = siteEllipsoids(parsed.sites);
        const fromStructure = buildElementColors(structure.elements, structure.elementCounts);
        const fromSites = buildElementColors(siteSpecies(parsed.sites), speciesCounts(sites));
        expect(fromSites).toEqual(fromStructure);
        expect(fromStructure).toEqual({ Ga: ELEMENT_COLORS.Ga, Se: ELEMENT_COLORS.Se, Ta: ELEMENT_COLORS.Ta });
    });
});
