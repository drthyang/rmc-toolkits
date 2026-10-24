// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// The holdout scan's rules, pinned on made-up cases (elements, lengths and
// directions that no validation case uses): what is a hit, and the false
// positives it must not raise. split.test.js runs it on the real test set.

import { describe, expect, it } from 'vitest';
import { describeHit, holdoutTokens, NEAR_WORDS, readText, scanText } from './holdoutScan';

const OFFCENTRE = { id: 'test-made-up', kind: 'offcentre', params: { elements: ['Eu', 'Lu', 'N'], a: 4.44, axis: [1, 1, 2], delta: 0.27 } };
const OFFCENTRE_DEV = { id: 'dev-made-up', kind: 'offcentre', params: { elements: ['La', 'Al', 'N'], a: 4.44, axis: [0, 0, 1], delta: 0.31 } };
const tokensOf = (item, devs) => holdoutTokens(item, devs);
const kindsFound = (text, tokens) => scanText(readText(text), tokens).map((hit) => `${hit.kind}:${hit.name}`);

describe('holdout tokens', () => {
    it('come from the parameters, and keep only what differs from the dev case of the same kind', () => {
        const tokens = tokensOf(OFFCENTRE, [OFFCENTRE_DEV]);
        expect(tokens.map((token) => `${token.kind}:${token.name}`).sort()).toEqual([
            'composition:elements', 'direction:axis', 'id:id', 'magnitude:delta'
        ]);
        // The same lattice parameter as the dev case is the kind's physics, not the test case.
        expect(tokens.some((token) => token.name === 'a')).toBe(false);
        // A dev case of another kind does not make an item ordinary.
        expect(tokensOf(OFFCENTRE, [{ ...OFFCENTRE_DEV, kind: 'thermal_null' }]).some((token) => token.name === 'a')).toBe(true);
    });

    it('a direction is compared as a cubic family, a wavevector up to sign, order and a reciprocal-lattice vector', () => {
        const zone = (axis, q) => ({ id: 'x', kind: 'zone_boundary', params: { elements: ['Eu', 'Lu', 'N'], a: 4, axis, q, delta: 0.2, decoy: { axis: [0, 1, 0], q: [0.125, 0, 0], amp: 0.1 } } });
        const tokens = tokensOf(zone([1, 0, 0], [0.75, 0, 0]), [zone([0, 0, -1], [0, 0.25, 0])]);
        expect(tokens.filter((token) => token.kind === 'direction' || token.kind === 'wavevector')).toEqual([]);
    });

    it('a structure word counts only when no dev case uses it and it is not a formula', () => {
        const sites = (structure) => ({ id: 'x', kind: 'sites_only', params: { structure, a: 6, sites: [{ el: 'Eu' }, { el: 'Lu' }] } });
        const words = tokensOf(sites('spinel EuLu2'), [sites('rocksalt LuN'), sites('spinel-like EuN')])
            .filter((token) => token.kind === 'structure').map((token) => token.word);
        expect(words).toEqual(['spinel']);
        expect(tokensOf(sites('rocksalt EuLu2'), [sites('Rocksalt LuN')]).filter((token) => token.kind === 'structure')).toEqual([]);
    });

    it('every kind of case has its reader', () => {
        expect(() => holdoutTokens({ id: 'x', kind: 'nope', params: {} }, [])).toThrow(/no tokens for case kind 'nope'/);
    });

    it('a judgment dev case (plan 1.6 ii) never removes a token: it is an option of the kind, not its physics', () => {
        const names = (tokens) => tokens.map((token) => `${token.kind}:${token.name}`).sort();
        // The test case's own delta, axis and lattice parameter, in a judgment dev case.
        const judged = { ...OFFCENTRE_DEV, id: 'dev-made-up-judged', judgment: true, params: { ...OFFCENTRE.params, elements: ['La', 'Al', 'N'] } };
        expect(names(tokensOf(OFFCENTRE, [judged]))).toEqual(['composition:elements', 'direction:axis', 'id:id', 'magnitude:a', 'magnitude:delta']);
        expect(tokensOf(OFFCENTRE, [OFFCENTRE_DEV, judged])).toEqual(tokensOf(OFFCENTRE, [OFFCENTRE_DEV]));
        // Nor does its structure word make a test case's word ordinary.
        const sites = (structure, extra = {}) => ({ id: 'x', kind: 'sites_only', ...extra, params: { structure, a: 6, sites: [{ el: 'Eu' }, { el: 'Lu' }] } });
        expect(tokensOf(sites('spinel EuLu2'), [sites('spinel LuN', { judgment: true })]).filter((token) => token.kind === 'structure')).toHaveLength(1);
    });

    it('the judgment subset\'s options and its own kind have their readers', () => {
        const names = (tokens) => tokens.filter((token) => token.kind !== 'id').map((token) => `${token.kind}:${token.name}`).sort();
        const weak = { id: 'test-made-up-weak', kind: 'zone_boundary', params: { elements: ['Eu', 'Lu', 'N'], a: 4.44, axis: [1, 1, 2], q: [0.25, 0.25, 0], delta: 0.017, decoy: null, weak: true } };
        expect(names(tokensOf(weak, []))).toEqual(['composition:elements', 'direction:axis', 'magnitude:a', 'magnitude:delta', 'wavevector:q']);
        const jt = { id: 'test-made-up-jt', kind: 'jahn_teller', params: { elements: ['Eu', 'Lu', 'N'], a: 9.1, short: 1.71, long: 2.07, offCentre: { axis: [1, 1, 2], delta: 0.13 } } };
        expect(names(tokensOf(jt, []))).toEqual(expect.arrayContaining(['direction:offCentre.axis', 'magnitude:offCentre.delta', 'magnitude:short']));
        const octahedra = { id: 'test-made-up-octahedra', kind: 'elongated_clusters', params: { element: 'Lu', anion: 'N', a: 9.1, edge: 2.71, elongation: 0.37 } };
        const tokens = tokensOf(octahedra, []);
        expect(names(tokens)).toEqual(['cluster:cluster', 'composition:elements', 'magnitude:a', 'magnitude:edge', 'magnitude:elongation', 'magnitude:long edge']);
        // The long edge is √(r0² + (r0 + elongation)²), r0 = edge/√2: 2.98 Å here.
        expect(kindsFound('Lu₆ octahedra, edges 2.71 and 2.98 Å around Lu', tokens).sort()).toEqual(['cluster:cluster', 'magnitude:edge', 'magnitude:long edge']);
    });
});

describe('the holdout scanner', () => {
    const tokens = tokensOf(OFFCENTRE, [OFFCENTRE_DEV]);

    it('finds the id anywhere, even inside a longer word', () => {
        expect(kindsFound('see test-made-up__B__r1.json', tokens)).toEqual(['id:id']);
    });

    it('a composition is a formula with exactly its elements, or a chain naming each once', () => {
        expect(kindsFound('EuLuN3', tokens)).toEqual(['composition:elements']);
        expect(kindsFound('Eu₂LuN₄ and (Eu2LuN4)', tokens)).toEqual(['composition:elements', 'composition:elements']);
        expect(kindsFound('the Eu–Lu–N angle; N-Lu-Eu', tokens)).toEqual(['composition:elements', 'composition:elements']);
        // Not hits: bare symbols, a subset, a superset, a chain that repeats an element.
        expect(kindsFound('Eu, Lu and N', tokens)).toEqual([]);
        expect(kindsFound('EuLu, LuN2 and EuLuN3O', tokens)).toEqual([]);
        expect(kindsFound('Lu–N–Lu and N–Lu–N–Eu–N', tokens)).toEqual([]);
        expect(kindsFound('Europium lutetium nitride', tokens)).toEqual([]);
    });

    it('a two-element composition is a formula only: its A–X–A chain is a bond angle', () => {
        const pair = tokensOf({ id: 'x', kind: 'thermal_null', params: { elements: ['Eu', 'N', 'N'], a: 4 } }, []);
        expect(kindsFound('EuN', pair)).toEqual(['composition:elements']);
        expect(kindsFound('Eu–N–Eu and N–Eu–N', pair)).toEqual([]);
    });

    it('a length counts within NEAR_WORDS words of its element, as written, with trailing zeros or not', () => {
        expect(NEAR_WORDS).toBe(6);
        expect(kindsFound('Lu is off-centred by ±0.27 Å', tokens)).toEqual(['magnitude:delta']);
        expect(kindsFound('a shift of 0.270 Å (Lu)', tokens)).toEqual(['magnitude:delta']);
        expect(kindsFound('in EuLuN3 the shift is 0.27 Å', tokens)).toEqual(['composition:elements', 'magnitude:delta']);
        expect(kindsFound('the Lu–N bond, 0.27 Å', tokens)).toEqual(['magnitude:delta']);
        expect(kindsFound('Lu one two three four five 0.27', tokens)).toEqual(['magnitude:delta']);
        // Not hits: too far, another element, other digits.
        expect(kindsFound('Lu one two three four five six 0.27', tokens)).toEqual([]);
        expect(kindsFound('La is off-centred by 0.27 Å', tokens)).toEqual([]);
        expect(kindsFound('Lu: 0.271, 10.27, 0.27.1 and 0.027', tokens)).toEqual([]);
        expect(kindsFound('a shift of 0.27 Å', tokens)).toEqual([]);
    });

    it('a direction counts as its family next to its element, however it is bracketed or spaced', () => {
        expect(kindsFound('Lu along [1 1 2]', tokens)).toEqual(['direction:axis']);
        expect(kindsFound('Lu along <211>', tokens)).toEqual(['direction:axis']);
        expect(kindsFound('Lu along ⟨1, −1, 2⟩', tokens)).toEqual(['direction:axis']);
        // A spaced vector is one word: "[1 1 2]" does not push the length out of reach.
        expect(kindsFound('Lu 0.10 Å [1 1 2] for 0.27', tokens)).toEqual(['magnitude:delta', 'direction:axis']);
        expect(kindsFound('Lu along [1 1 0] and [0 0 1]', tokens)).toEqual([]);
        expect(kindsFound('a vector [1 1 2] with no element', tokens)).toEqual([]);
    });

    it('a wavevector counts up to sign, order and G, next to its element', () => {
        const zone = tokensOf({
            id: 'x', kind: 'zone_boundary',
            params: { elements: ['Eu', 'Lu', 'N'], a: 4, axis: [0, 0, 1], q: [0.25, 0.25, 0], delta: 0.2, decoy: { axis: [0, 1, 0], q: [0.125, 0, 0], amp: 0.1 } }
        }, [{ id: 'y', kind: 'zone_boundary', params: { elements: ['La', 'Al', 'N'], a: 4, axis: [0, 0, 1], q: [0, 0, 0.5], delta: 0.2, decoy: { axis: [0, 1, 0], q: [0.125, 0, 0], amp: 0.1 } } }]);
        const found = (text) => kindsFound(text, zone).filter((kind) => kind.startsWith('wavevector'));
        expect(found('Lu at q = (¼, ¼, 0)')).toEqual(['wavevector:q']);
        expect(found('Lu at q = (0, 1/4, -1/4)')).toEqual(['wavevector:q']);
        expect(found('Lu at q = (0.75, 0.25, 0)')).toEqual(['wavevector:q']);
        expect(found('Lu at q = (½, 0, 0)')).toEqual([]);
        expect(found('q = (¼, ¼, 0) with no element near it at all here')).toEqual([]);
    });

    it('a cluster is its own formula, not part of a larger one', () => {
        const dimer = tokensOf({ id: 'x', kind: 'rigid_dimer', params: { element: 'Lu', a: 5, axis: 0, bond: 2.4, sigmaCommon: 0.1 } }, []);
        const found = (text) => kindsFound(text, dimer).filter((kind) => kind.startsWith('cluster'));
        expect(found('Lu₂ dimers')).toEqual(['cluster:cluster']);
        expect(found('rigid Lu2-pairs')).toEqual(['cluster:cluster']);
        expect(found('Lu2O3, Lu, Lu3 and EuLu2')).toEqual([]);
    });

    it('a structure word is a whole word, in any case', () => {
        const sites = tokensOf({ id: 'x', kind: 'sites_only', params: { structure: 'spinel EuLu2', a: 6, sites: [{ el: 'Eu' }, { el: 'Lu' }] } }, []);
        const found = (text) => kindsFound(text, sites).filter((kind) => kind.startsWith('structure'));
        expect(found('a Spinel, the spinel')).toEqual(['structure:structure', 'structure:structure']);
        expect(found('spinels and inverse-spinelike')).toEqual([]);
    });

    it('a seed counts next to the word "seed", as a whole number, wherever its case', () => {
        const seeded = tokensOf({ ...OFFCENTRE, params: { ...OFFCENTRE.params, seed: 987 } }, [{ ...OFFCENTRE_DEV, params: { ...OFFCENTRE_DEV.params, seed: 987 } }]);
        // A seed is never the kind's physics: the dev case's equal seed does not remove it.
        expect(seeded.filter((token) => token.kind === 'seed')).toEqual([{ case: 'test-made-up', kind: 'seed', name: 'seed', value: 987 }]);
        const found = (text) => kindsFound(text, seeded).filter((kind) => kind.startsWith('seed'));
        expect(found('seed 987')).toEqual(['seed:seed']);
        expect(found('seeds 985, 986, 987 in order')).toEqual(['seed:seed']);
        expect(found('{ "seed": 987 }')).toEqual(['seed:seed']);
        expect(found('the Seed was 987.')).toEqual(['seed:seed']);
        // Not hits: no "seed" near it, or another number.
        expect(found('987 atoms')).toEqual([]);
        expect(found('seed one two three four five six 987')).toEqual([]);
        expect(found('seed 9870, 1987, 2.987 and 987.5')).toEqual([]);
        expect(found('seedling 987')).toEqual([]);
        // The digits that end an identifier are not the number (TEST_SET_SHA256 beside "seed").
        expect(found('each seed is pinned by TEST_SET_SHA987')).toEqual([]);
    });

    it('a supercell counts next to its element, written as N³, N × N × N or [N, N, N]; the dev twin\'s does not', () => {
        const cubic = (supercell) => ({ ...OFFCENTRE, params: { ...OFFCENTRE.params, supercell } });
        const tokens = tokensOf(cubic([9, 9, 9]), [{ ...OFFCENTRE_DEV, params: { ...OFFCENTRE_DEV.params, supercell: [5, 5, 5] } }]);
        const found = (text) => kindsFound(text, tokens).filter((kind) => kind.startsWith('supercell'));
        expect(found('Lu in a 9³ box')).toEqual(['supercell:supercell']);
        expect(found('Lu, 9 × 9 × 9 cells')).toEqual(['supercell:supercell']);
        expect(found('Lu supercell 9x9x9')).toEqual(['supercell:supercell']);
        expect(found('Lu, supercell: [9, 9, 9]')).toEqual(['supercell:supercell']);
        // Not hits: far from its element, another size, or only next to the word.
        expect(found('a 9³ supercell')).toEqual([]);
        expect(found('Lu in a 19³ box, 9 × 9 × 8 and 9.5³')).toEqual([]);
        expect(tokensOf(cubic([5, 5, 5]), [{ ...OFFCENTRE_DEV, params: { ...OFFCENTRE_DEV.params, supercell: [5, 5, 5] } }])
            .some((token) => token.kind === 'supercell')).toBe(false);
    });

    it('a thermal width the case sets, and a breathing sd as a fraction or in per cent, count as magnitudes next to the element', () => {
        const breathing = { id: 'test-made-up-breathing', kind: 'breathing', params: { shape: 'triangle', element: 'Lu', anion: 'N', a: 7.7, edge: 2.91, breathe: 0.027, sigma: 0.043 } };
        const dev = { id: 'dev-made-up-breathing', kind: 'breathing', params: { shape: 'triangle', element: 'La', anion: 'N', a: 7.1, edge: 2.71, breathe: 0.045, sigma: 0.043 } };
        const tokens = tokensOf(breathing, [dev]);
        const names = tokens.filter((token) => token.kind === 'magnitude').map((token) => token.name).sort();
        // The dev twin's thermal width is the kind's: not a token; the breathing sd is.
        expect(names).toEqual(['a', 'breathe', 'breathe %', 'edge']);
        const found = (text) => kindsFound(text, tokens).filter((kind) => kind.startsWith('magnitude'));
        expect(found('Lu₃ triangles breathe with sd 0.027')).toEqual(['magnitude:breathe']);
        expect(found('Lu₃ triangles, sd 2.7 %')).toEqual(['magnitude:breathe %']);
        expect(names.includes('sigma')).toBe(false);
        const jt = { id: 'test-made-up-jt-thermal', kind: 'jahn_teller', params: { elements: ['Eu', 'Lu', 'N'], a: 9.1, short: 1.71, long: 2.07, sigma: 0.047 } };
        expect(kindsFound('Lu–N thermal sd 0.047 Å', tokensOf(jt, []))).toEqual(['magnitude:sigma']);
    });

    it('reports a hit by file, line, case and kind with what matched masked', () => {
        const text = 'line one\nLu is off-centred by 0.27 Å in EuLuN3\n';
        const read = readText(text);
        const hits = scanText(read, tokens);
        const lines = hits.map((hit) => describeHit(read, hit, hits, 'docs/x.md'));
        expect(lines).toEqual([
            'docs/x.md:2 test-made-up composition (elements): ■ is off-centred by ■ Å in ■',
            'docs/x.md:2 test-made-up magnitude (delta): ■ is off-centred by ■ Å in ■'
        ]);
        lines.forEach((line) => {
            expect(line).not.toContain('0.27');
            expect(line).not.toContain('EuLuN3');
        });
    });
});
