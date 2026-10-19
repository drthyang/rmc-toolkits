// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// The verdict's words as claims (verdict.js VERDICT_PHRASINGS, adapter 2 of
// the verification's `phrasings` seam). The app's own claim catalogue reads
// prose only, since the verification module (behaviour change: before, it
// read the schema's words itself); the harness passes these to every turn
// (`config.phrasings`), so a VERDICT line still needs its analysis as the
// same claim in prose would. Each phrasing is matched apart from the app's,
// and a claim holds when either does: the readings the harness's tests pin
// (harness.test.js, unchanged) are the union of the two.

import { describe, expect, it } from 'vitest';
import { CLAIMS, createVerification } from '../../src/llm/copilot/verification';
import { FEATURE_KINDS, VERDICT_PHRASINGS } from '../verdict';

// The claims an answer makes: judged with nothing run, each is an evidence finding.
const claims = (text, config = {}) => createVerification(config).judge(text).done.audit
    .filter((finding) => finding.rule === 'evidence').map((finding) => finding.id);
const harness = (text) => claims(text, { phrasings: VERDICT_PHRASINGS });
// The claims whose verdict-word phrasing matches the text.
const adapterMatches = (text) => Object.entries(VERDICT_PHRASINGS)
    .filter(([, { asserted = [], literal = [] }]) => [...asserted, ...literal].some((pattern) => text.search(pattern) >= 0))
    .map(([id]) => id);
// In CLAIMS order, as the gate reads them.
const union = (...lists) => CLAIMS.map((claim) => claim.id).filter((id) => lists.some((list) => list.includes(id)));

const verdict = (features, rest = '') => `VERDICT {"features": ${JSON.stringify(features)}${rest}}`;
// The kind each feature of the verdict is read as: the claim that needs its analysis.
const KIND_CLAIMS = {
    off_centre: 'split-site',
    modulation: 'order-between-cells',
    cluster_distortion: 'clusters',
    cluster_breathing: 'clusters',
    rigid_cluster: 'clusters',
    jahn_teller: 'coordination',
    symmetry_lowering: 'lowering'
};

describe('the verdict\'s words, read for the evidence gate (adapter 2)', () => {
    it('maps every feature kind of the schema to the claim that needs its analysis', () => {
        expect(Object.keys(KIND_CLAIMS).sort()).toEqual([...FEATURE_KINDS].sort());
        const ids = CLAIMS.map((claim) => claim.id);
        Object.values(KIND_CLAIMS).forEach((id) => expect(ids).toContain(id));
        Object.keys(VERDICT_PHRASINGS).forEach((id) => expect(ids).toContain(id));
        // The phrasings the turn accepts: a known claim, RegExp patterns, the asserted ones global.
        expect(() => createVerification({ phrasings: VERDICT_PHRASINGS })).not.toThrow();
    });

    it.each(FEATURE_KINDS)('a %s feature needs its analysis in the harness; only "modulation" is a word the app reads itself', (kind) => {
        const text = verdict([{ kind, element: 'Ti', direction_uvw: [0, 0, 1], magnitude_A: null, q: null }]);
        expect(harness(text)).toEqual([KIND_CLAIMS[kind]]);
        // The negative case: the app's own catalogue reads no schema word.
        expect(claims(text)).toEqual(kind === 'modulation' ? ['order-between-cells'] : []);
    });

    it('an empty feature list is "no order" in the harness, anywhere in the text, and nothing to the app', () => {
        const text = verdict([], ', "space_group": "Pm-3m"');
        expect(harness(text)).toEqual(['no-order']);
        expect(claims(text)).toEqual([]);
        // The VERDICT-word cases of the gate's tests, as they were pinned.
        expect(harness('VERDICT {"features": [{"kind": "off_centre", "element": "Nb"}]}')).toEqual(['split-site']);
        const moClusters = { elements: ['Mo', 'O'], shell: (centre, partner) => (centre === partner ? { cn: centre === 'Mo' ? 2 : 6 } : null) };
        expect(claims('VERDICT {"features": [], "space_group": "Pm-3m"}', { catalog: moClusters, phrasings: VERDICT_PHRASINGS })).toEqual(['no-order']);
    });

    // The final answers of the harness's pinned runs (harness.test.js,
    // fitCase.test.js): what the harness's turn reads in each is the app's
    // reading or the verdict words' — each matched apart.
    const PINNED = [
        ['dev-jahn-teller: coordination in prose, "no order" only in the verdict',
            'Cu is four-coordinated (square planar), with the next F shell at 2.25 Å.\n\nVERDICT {"features": [], "space_group": "Pm-3m", "chiral": false, "centrosymmetric": true, "polar": false}',
            ['coordination'], ['coordination', 'no-order']],
        ['dev-zone-boundary: a split site in prose and in the verdict',
            'Ti shows random two-site disorder along c (a split site).\n\nVERDICT {"features": [{"kind": "off_centre", "element": "Ti", "direction_uvw": [0, 0, 1], "magnitude_A": 0.1}], "space_group": "Pm-3m", "chiral": false, "centrosymmetric": true, "polar": false}',
            ['split-site'], ['split-site']],
        ['dev-chiral: thermal motion only, and an empty verdict',
            'Thermal motion only.\n\nVERDICT {"features": [], "space_group": "P2_13", "chiral": true, "centrosymmetric": false, "polar": false}',
            ['no-order'], ['no-order']],
        ['dev-offcentre: "No order." and an empty verdict',
            'No order.\n\nVERDICT {"features": [], "space_group": "Pm-3m", "chiral": false, "centrosymmetric": true, "polar": false}',
            ['no-order'], ['no-order']],
        ['dev-breathing: rigid units in prose and in the verdict',
            'The Mo3 triangles move as rigid units.\n\nVERDICT {"features": [{"kind": "rigid_cluster", "element": "Mo"}], "space_group": null, "chiral": null, "centrosymmetric": null, "polar": null}',
            ['clusters'], ['clusters']],
        ['dev-symmetry-lowering: a quantified shift in prose and in the verdict',
            'Ta is shifted 0.061 Å along [0, 0, 1]: the average structure is P4mm, not Pm-3m.\n\nOutcome: the shift the card hides.\n\n'
                + 'VERDICT {"features": [{"kind": "symmetry_lowering", "element": "Ta", "direction_uvw": [0, 0, 1], "magnitude_A": 0.061, "q": null}], "space_group": "P4mm", "chiral": false, "centrosymmetric": false, "polar": true}',
            ['lowering'], ['lowering']],
        ['dev-fit-unconverged: "no order beyond the average structure" and an empty verdict',
            'The run has converged and the fit is good. Thermal motion only: no order beyond the average structure.\n\nVERDICT {"features": [], "space_group": "Fm-3m", "chiral": false, "centrosymmetric": true, "polar": false}',
            ['no-order'], ['no-order']]
    ];

    it.each(PINNED)('%s', (_, text, app, pinned) => {
        expect(claims(text)).toEqual(app);
        expect(harness(text)).toEqual(pinned);
        // Matched apart and OR-ed: the harness's reading is the union of the app's and the verdict words'.
        expect(harness(text)).toEqual(union(claims(text), adapterMatches(text)));
    });
});
