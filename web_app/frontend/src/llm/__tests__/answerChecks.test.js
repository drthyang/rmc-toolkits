// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// The checks on the answer (the verification's facts catalogue,
// copilot/verification/facts.js) at the interface: each tool's facts
// recorded as the step of the tool that returns them, the run context given,
// then the answer judged. No claim row is offered (offered: []), so the
// findings are the checks' alone. Each case keeps its text and its expected
// findings from when the checks were answerAudit.js; a context given then as
// the checks' own facts is given here as the run context that yields them.

import { describe, expect, it } from 'vitest';
import { toolResult } from '../copilot/tools';
import { createVerification, FACTS, normalizeSymbol } from '../copilot/verification';
import { pointGroupFacts } from '../context/pointGroups';

// The tool that returns a facts object's key (FACTS `from`).
const producer = (facts) => FACTS.find((row) => row.key in facts)?.from[0] ?? 'some_tool';
const audit = (text, { facts = [], context = null } = {}) => {
    const verification = createVerification({ context, offered: [] });
    facts.forEach((item) => verification.record({ name: producer(item), ...toolResult({ facts: item }) }, () => ({ data: null })));
    return verification.judge(text).done.audit;
};
const rules = (text, options) => audit(text, options).map((finding) => finding.rule);

describe('the convergence check', () => {
    const improving = { facts: [{ convergence: 'improving' }] };

    it('flags "converged" against an improving run, from a tool or from the context', () => {
        expect(rules('The run has converged: the χ² history has reached a plateau.', improving)).toEqual(['convergence']);
        expect(rules('Yes — it is converged.', { context: { convergence: { watchdog_status: 'diverging' } } })).toEqual(['convergence']);
        expect(audit('It looks converged.', improving)[0].text).toMatch(/still improving .* it has not converged/);
    });

    it('stays quiet on a negated claim, a question, an acknowledged trend, or a converged run', () => {
        for (const text of [
            'The run has not converged; it is still improving.',
            'It hasn\'t converged yet.',
            'Has the run converged? No: ln χ² still falls.',
            'Although the curve looks converged, the watchdog says improving.',
            'The run is far from converged.'.replace('far from', 'not'),
            'Convergence would need a flat history.'
        ]) {
            expect(rules(text, improving)).toEqual([]);
        }
        expect(rules('The run has converged.', { facts: [{ convergence: 'converged' }] })).toEqual([]);
        expect(rules('The run has converged.', {})).toEqual([]);
    });
});

describe('claims, not mentions (review 2026-10-12 #3 and #15)', () => {
    const improving = { facts: [{ convergence: 'improving' }] };
    const none = { facts: [{ direction: { subject: 'site 23', significant: false, peak_sigma: 1.2, map_sigma: 0.8 } }] };

    it('does not read a hypothetical, conditional or a question as a claim', () => {
        for (const text of [
            'To check whether the run has converged, I read the χ² history: it is still falling.',
            'Keep refining until χ² plateaus.',
            'If the run had converged, the χ² history would be flat.',
            'Once it has converged, compare the g(r) again.',
            'A converged run would show a flat history; this one falls.'
        ]) {
            expect(rules(text, improving), text).toEqual([]);
        }
        expect(rules('I tested whether site 23 has a preferred direction: the peak is 1.2σ.', none)).toEqual([]);
    });

    it('flags a claim that follows a negation in another part of the sentence', () => {
        expect(rules('With no remaining drift in χ², the run has converged.', improving)).toEqual(['convergence']);
        expect(rules('There is no drift left, and the run has converged.', improving)).toEqual(['convergence']);
        expect(rules('The run has converged; χ² stopped improving.', improving)).toEqual(['convergence']);
        expect(rules('Without any doubt, site 23 has a preferred displacement direction.', none)).toEqual(['direction']);
    });
});

describe('the symmetry checks', () => {
    // F-43m above P-42₁m, P2₁3 and Pm-3m on the ladder: F-43m is the parent
    // (the average group and the ladder's top), as the four groups given
    // without a parent made it before.
    const context = {
        symmetry: {
            space_group: 'F-43m', point_group: '-43m',
            ladder: [
                { sg: 'F-43m', pg: '-43m', n_ops: 96 }, { sg: 'P-42_1m', pg: '-42m', n_ops: 8 },
                { sg: 'P2_13', pg: '23', n_ops: 12 }, { sg: 'Pm-3m', pg: 'm-3m', n_ops: 48 }
            ]
        }
    };

    it('reads a symbol however it is written', () => {
        expect(normalizeSymbol('P-42₁m')).toBe('p-421m');
        expect(normalizeSymbol('P−4 2_1 m')).toBe('p-421m');
        expect(normalizeSymbol('$P\\bar{4}2_1m$')).toBe('p-421m');
        expect(normalizeSymbol('F-43m (No. 216)')).toBe('f-43m');
    });

    it('flags the spot check\'s error: P-42₁m called chiral, F-43m said to lose its centrosymmetry', () => {
        const text = 'At 0.02 Å the group is P-42₁m, which is chiral. Going from F-43m to P-42₁m breaks the centrosymmetry.';
        const found = audit(text, { context });
        expect(found.map((finding) => finding.rule)).toEqual(['chirality', 'inversion-lost']);
        expect(found[0].text).toMatch(/P-42_1m \(-42m\) chiral/);
        expect(found[1].text).toMatch(/none of F-43m \(-43m\), P-42_1m \(-42m\) has one/);
        expect(rules('F-43m is centrosymmetric.', { context })).toEqual(['centrosymmetry']);
        expect(rules('The centre of symmetry of $P\\bar{4}2_1m$ is lost at low tolerance.', { context })).toEqual(['inversion-lost']);
    });

    it('stays quiet on correct statements, negations and Sohncke or centrosymmetric groups', () => {
        for (const text of [
            'P-42₁m is achiral, and neither F-43m nor P-42₁m is centrosymmetric.',
            'F-43m is non-centrosymmetric; P-42₁m is not chiral either.',
            'P2₁3 is a Sohncke group, so the structure can be chiral.',
            'P2₁3 is chiral.',
            'Pm-3m is centrosymmetric, and lowering it to F-43m breaks the inversion.',
            'Does P-42₁m lose the centrosymmetry?'
        ]) {
            expect(rules(text, { context })).toEqual([]);
        }
    });

    it('accepts a lost inversion when the parent has one, even if a tool found a group without (review 2026-10-12 #4)', () => {
        // The average group Pm-3m has inversion; a tool found P4mm, which has
        // none. "Loses its inversion centre" is right, and no group is named.
        const facts = [{ groups: [{ symbol: 'P4mm', ...pointGroupFacts('4mm') }] }];
        const lowered = { facts, context: { symmetry: { space_group: 'Pm-3m', point_group: 'm-3m' } } };
        expect(rules('At 0.01 Å the local structure loses its inversion centre.', lowered)).toEqual([]);
        expect(rules('Going to P4mm breaks the inversion symmetry of the average structure.', lowered)).toEqual([]);
        // With a parent that has no inversion, the claim is still caught.
        expect(rules('At 0.01 Å the local structure loses its inversion centre.', {
            facts, context: { symmetry: { space_group: 'F-43m', point_group: '-43m' } }
        })).toEqual(['inversion-lost']);
    });

    it('does not blame the parent for a subgroup it does not know (review 2026-10-12 #14)', () => {
        const parentOnly = { context: { symmetry: { space_group: 'P-42_1m', point_group: '-42m' } } };
        expect(rules('The parent group is P-42₁m; its subgroup P2₁2₁2 is chiral.', parentOnly)).toEqual([]);
        // A clause that names no group is still read with its sentence's.
        expect(rules('The group is P-42₁m; it is chiral.', parentOnly)).toEqual(['chirality']);
    });

    it('takes the groups from a tool\'s facts too', () => {
        const facts = [{ groups: [{ symbol: 'Cmm2', ...pointGroupFacts('mm2') }] }];
        expect(rules('Cmm2 is centrosymmetric.', { facts })).toEqual(['centrosymmetry']);
        expect(rules('Cmm2 is polar.', { facts })).toEqual([]);
    });

    it('reads the parents from the run context: the average group and the ladder\'s top, or its two ends without counts', () => {
        // The watchdog's status and every rung's group, the space group as written without its number.
        const ladder = {
            convergence: { watchdog_status: 'improving' },
            symmetry: { space_group: 'F-43m (No. 216)', point_group: '-43m', ladder: [{ sg: 'P-42_1m', pg: '-42m' }, { sg: 'P1', pg: '1' }] }
        };
        expect(rules('The run has converged.', { context: ladder })).toEqual(['convergence']);
        expect(rules('P1 is chiral.', { context: ladder })).toEqual([]);
        expect(audit('The inversion centre is lost.', { context: ladder })[0].text)
            .toBe('The answer says an inversion centre is lost, but none of F-43m (-43m), P-42_1m (-42m), P1 (1) has one: there is no inversion to lose.');
        // With operation counts, the ladder's top is the rung with the most …
        const counted = { symmetry: { ladder: [{ sg: 'P4mm', pg: '4mm', n_ops: 8 }, { sg: 'Pm-3m', pg: 'm-3m', n_ops: 48 }, { sg: 'P1', pg: '1', n_ops: 1 }] } };
        expect(rules('The inversion centre is lost.', { context: counted })).toEqual([]);
        // … without them, the ladder's two ends.
        const uncounted = { symmetry: { ladder: [{ sg: 'P4mm', pg: '4mm' }, { sg: 'Pm-3m', pg: 'm-3m' }, { sg: 'P1', pg: '1' }] } };
        expect(audit('The inversion centre is lost.', { context: uncounted })[0].text).toMatch(/none of P4mm \(4mm\), P1 \(1\) has one/);
        expect(rules('The run has converged.', { context: null })).toEqual([]);
    });
});

describe('the preferred-direction check', () => {
    const none = { facts: [{ direction: { subject: 'Ga', significant: false, peak_sigma: -0.887, map_sigma: 0.561 } }] };

    it('flags a preferred direction no analysis of the question found', () => {
        const found = audit('The Ga atoms have a preferred displacement direction along [111].', none);
        expect(found.map((finding) => finding.rule)).toEqual(['direction']);
        expect(found[0].text).toMatch(/\(Ga: peak -0.887σ, map 0.561σ\)/);
    });

    it('stays quiet on "no preferred direction", or when an analysis found one', () => {
        expect(rules('The Ga atoms do not exhibit a significant preferred displacement direction.', none)).toEqual([]);
        expect(rules('There is no preferred direction.', none)).toEqual([]);
        expect(rules('Site 1 has a preferred displacement direction along [111].', {
            facts: [...none.facts, { direction: { subject: 'site 1', significant: true, peak_sigma: 4.2, map_sigma: 3.1 } }]
        })).toEqual([]);
        expect(rules('A preferred direction along [111].', {})).toEqual([]);
    });
});

describe('the polar check (local-symmetry skill) and the superstructure check (validation 2026-10-14)', () => {
    const cubic = { facts: [{ groups: [{ symbol: 'Pm-3m', ...pointGroupFacts('m-3m') }] }] };
    const tetragonal = { facts: [{ groups: [{ symbol: 'P4mm', ...pointGroupFacts('4mm') }] }] };

    it('flags "polar" for a group whose point group is not, and nothing else', () => {
        expect(rules('The average structure, Pm-3m, is polar along c.', cubic)).toEqual(['polar']);
        for (const text of [
            'Pm-3m is non-polar and centrosymmetric.',
            'Pm-3m is not polar.',
            'Is Pm-3m polar? No.',
            'Pm-3m has no polar axis.'
        ]) {
            expect(rules(text, cubic)).toEqual([]);
        }
        expect(rules('P4mm is polar along c.', tetragonal)).toEqual([]);
    });

    it('a JSON key is a field name, not a claim (a structured verdict under the answer)', () => {
        const json = 'VERDICT {"features": [], "space_group": "Pm-3m", "chiral": false, "centrosymmetric": true, "polar": false}';
        expect(rules(json, cubic)).toEqual([]);
        expect(rules(`Pm-3m is polar.\n\n${json}`, cubic)).toEqual(['polar']);
    });

    it('flags a superstructure claim only when every modulation analysis found none', () => {
        const none = { facts: [{ modulation: { elements: ['Ti'], ordered: [] } }] };
        const found = { facts: [{ modulation: { elements: ['Ti'], ordered: [{ element: 'Ti', q: [0, 0, 0.5] }] } }] };
        expect(rules('Ti forms a superstructure that doubles the cell along c.', none)).toEqual(['superstructure']);
        expect(rules('The Ti displacements show zone-boundary order at (0, 0, ½).', none)).toEqual(['superstructure']);
        for (const text of [
            'There is no superstructure: no zone-boundary point stands out.',
            'I tested whether a superstructure doubles the cell.',
            'Without a superstructure, the split is random two-site disorder.'
        ]) {
            expect(rules(text, none)).toEqual([]);
        }
        expect(rules('Ti forms a superstructure that doubles the cell along c.', found)).toEqual([]);
        expect(rules('Ti forms a superstructure.', {})).toEqual([]);
    });
});

describe('cluster and coordination claims (review L5)', () => {
    const breathing = { facts: [{ cluster: { element: 'Mo', members: 3, distorted: false, static: false, breathing: true, rigid: false, split: false } }] };
    const jt = { facts: [{ cluster: { element: 'F', members: 6, distorted: true, static: false, breathing: false, rigid: false, split: true } }] };

    it('flags a cluster distortion, breathing or rigidity every analysis of that element contradicts', () => {
        expect(rules('The Mo3 triangles are distorted beyond thermal motion.', breathing)).toEqual(['cluster']);
        expect(rules('The Mo clusters move as rigid units.', breathing)).toEqual(['cluster']);
        expect(audit('The Mo3 triangles are distorted.', breathing)[0].text).toMatch(/^The answer says the Mo3 clusters are distorted, but analyze_clusters found no shape distortion/);
        const rigidDimers = { facts: [{ cluster: { element: 'V', members: 2, distorted: false, static: false, breathing: false, rigid: true, split: false } }] };
        expect(rules('The V2 dimers breathe.', rigidDimers)).toEqual(['cluster']);
    });

    it('stays quiet on what an analysis found, a negation, split sites, and a polyhedron of another element', () => {
        for (const text of [
            'The Mo3 clusters breathe: their size fluctuates collectively.',
            'The Mo3 triangles are not distorted; their size fluctuates.',
            'The TaSe6 octahedra are distorted.',
            'I checked whether the Mo3 clusters are distorted.'
        ]) {
            expect(rules(text, breathing)).toEqual([]);
        }
        expect(rules('The F6 octahedra are distorted (Jahn–Teller).', jt)).toEqual([]);
        expect(rules('The Mo3 clusters are distorted.', {})).toEqual([]);
    });

    const shell = { facts: [{ coordination: [
        { centre: 'Cu', partner: 'F', window: [1.7, 2.52], cn: 6, whole_shell: true, confirmed: false },
        { centre: 'F', partner: 'Cu', window: [1.7, 2.52], cn: 1, whole_shell: true, confirmed: false }
    ] }] };

    it('flags a coordination number the whole first shell contradicts', () => {
        expect(rules('Cu is four-coordinated (square planar).', shell)).toEqual(['coordination']);
        expect(rules('Each Cu has CN 4.', shell)).toEqual(['coordination']);
        expect(rules('Cu sits in a square-planar coordination.', shell)).toEqual(['coordination']);
        expect(audit('Cu is 4-coordinated.', shell)[0].text).toMatch(/^The answer gives Cu a coordination of 4, but its whole first shell holds 6 F per Cu in 1\.7–2\.52 Å/);
    });

    it('stays quiet on the right count, a partial shell, a hypothetical, and a window the user chose', () => {
        for (const text of [
            'Cu is six-coordinated by F.',
            'Cu has an octahedral coordination with 4 short and 2 long bonds.',
            'Within the first peak, Cu is four-coordinated.',
            'In 1.7–2.1 Å Cu is four-coordinated.',
            'Cu would be four-coordinated if only the short bonds counted.',
            'The coordination number of Cu is not 4.'
        ]) {
            expect(rules(text, shell)).toEqual([]);
        }
        const chosen = { facts: [...shell.facts, { coordination: { centre: 'Cu', partner: 'F', window: [1.7, 2.1], cn: 4, whole_shell: false, confirmed: true } }] };
        expect(rules('Cu is four-coordinated.', chosen)).toEqual([]);
        expect(rules('Cu is four-coordinated.', {})).toEqual([]);
    });

    it('a general wavevector\'s order allows "superstructure", never "cell doubling"', () => {
        const general = { facts: [{ modulation: { elements: ['Ti'], ordered: [], general: [{ element: 'Ti', q: [0.375, 0, 0] }] } }] };
        expect(rules('Ti shows an incommensurate superstructure.', general)).toEqual([]);
        expect(rules('Ti order doubles the cell along a.', general)).toEqual(['superstructure']);
    });
});

describe('the long-wave check: order at a box-scale wavevector (the user\'s rule, 2026-10-24)', () => {
    // analyze_modulations' facts: the engine's bound (every |q_i| ≤ 0.25 is a
    // long wave), Ta's strongest mode a long wave, an ordered zone-boundary point.
    const facts = [{
        modulation: {
            elements: ['Ta'], ordered: [{ element: 'Ta', q: [0, 0, 0.5] }], general: [], long_wave_q: 0.25,
            long_waves: [{ element: 'Ta', q: [0.1, 0, 0] }]
        }
    }];

    it('flags order claimed at a long wave, by its q or in words', () => {
        for (const text of [
            'There is long-range order at q = (0.1, 0, 0): the strongest Ta mode, 6× the median.',
            'The Ta displacements are modulated with q = (0.1 0.1 0).',
            'Ta shows ordered displacements at q=[1/8, 0, 0].',
            'A superstructure with wavevector (0.25, 0, 0) appears in the Ta modes.',
            'Ta shows a long-wavelength modulation, which indicates incommensurate order.',
            'The box-scale wave is evidence of ordered Ta displacements.'
        ]) {
            expect({ text, rules: rules(text, { facts }) }).toEqual({ text, rules: ['long-wave'] });
        }
        expect(audit('Ta is ordered at q = (0.1, 0, 0).', { facts })[0].text).toBe('The answer claims order at q = (0.1, 0, 0), a box-scale '
            + 'wavevector (every |q_i| ≤ 0.25: a long wave), but a long wave is a correlation the finite box can\'t represent, so it says '
            + 'nothing about order at its own q; a point or mode that passes its own test at another q is judged on that test alone.');
    });

    it('stays quiet on the user\'s own wording, a denial, a long wave named without order, and order that passed its own test elsewhere', () => {
        for (const text of [
            'The strongest Ta mode is a long wave at q = (0.1, 0, 0): a correlation the finite box can\'t represent, so it says nothing about order at its own q.',
            'Ta\'s strongest mode is a long-wavelength wave, q = (0.1, 0, 0) at 6× the median: box-scale displacement correlations, not a crystallographic superstructure.',
            'The modulation at q = (0.1, 0, 0) is a long wave, not order.',
            'No order: the long waves are box-scale correlations.',
            'Ordered Ta displacements at q = (0, 0, ½): a superstructure doubling the cell along c.',
            'Order at q = (0.3, 0, 0) passed its own test.',
            'Ta displacements are ordered along [0 0 1].',
            'In order to judge the long waves, I ran analyze_modulations.'
        ]) {
            expect({ text, rules: rules(text, { facts }) }).toEqual({ text, rules: [] });
        }
    });

    it('needs an analyze_modulations of the question (its bound), as the superstructure check', () => {
        expect(rules('There is order at q = (0.1, 0, 0).', {})).toEqual([]);
        // The superstructure rule still judges a claim without an ordered point.
        expect(rules('Ta is ordered at q = (0.1, 0, 0) — a superstructure.', { facts: [{ modulation: { ...facts[0].modulation, ordered: [] } }] }))
            .toEqual(['superstructure', 'long-wave']);
    });
});

describe('the single-configuration check: a distribution quoted without saying it is one configuration\'s (the user\'s rule, 2026-10-24)', () => {
    const angles = [{ distribution: { kind: 'bond_angles', subject: 'the Se–Ga–Se bond-angle distribution', configurations: 1 } }];

    it('flags a distribution\'s numbers quoted without saying they come from one configuration', () => {
        for (const text of [
            'Se–Ga–Se angles: mean 109.4°, sd 4.1°, mode 109.5° with share 0.94.',
            'An octahedral split: 80 % cis near 89°, 20 % trans near 159°.',
            'The sin θ-corrected histogram peaks at 100–110°.',
            '3 948 of 4 000 Ga atoms have 4 neighbours, mean 3.99.',
            'The angle distribution is centred at 109.4° with a spread of 4.1°.'
        ]) {
            expect({ text, rules: rules(text, { facts: angles }) }).toEqual({ text, rules: ['single-configuration'] });
        }
        expect(audit('Se–Ga–Se angles: mean 109.4°, sd 4.1°.', { facts: angles })[0].text).toBe('The answer quotes the numbers of the Se–Ga–Se '
            + 'bond-angle distribution without saying they come from one configuration: every analysis here reads the run\'s one '
            + 'configuration, so its distribution carries no spread between runs; say so.');
    });

    it('stays quiet when the answer says one configuration, quotes no distribution, or no distribution analysis ran', () => {
        for (const text of [
            'In this one configuration the Se–Ga–Se angles have mean 109.4°, sd 4.1°.',
            'Se–Ga–Se angles: mean 109.4°, sd 4.1° (a single RMC configuration; independent runs would give the spread).',
            'From one configuration: 80 % cis near 89°, 20 % trans near 159°.',
            'The GaSe4 tetrahedron is regular.',
            'There are no near-180° angles.',
            'Ga–Se 2.434 Å, sd 0.095 Å.',
            'No preferred direction (peak −0.89σ, map 0.56σ).'
        ]) {
            expect({ text, rules: rules(text, { facts: angles }) }).toEqual({ text, rules: [] });
        }
        expect(rules('Se–Ga–Se angles: mean 109.4°, sd 4.1°.', {})).toEqual([]);
        expect(rules('Se–Ga–Se angles: mean 109.4°, sd 4.1°.', { facts: [{ coordination: { centre: 'Ga', partner: 'Se', window: [2.18, 2.79], cn: 4, whole_shell: true, confirmed: false } }] })).toEqual([]);
    });
});

describe('the order of the checks', () => {
    it('runs in FACTS order, and a lowering\'s shifts have no check yet', () => {
        expect(FACTS.map((row) => row.key)).toEqual(['convergence', 'groups', 'direction', 'modulation', 'cluster', 'coordination', 'distribution', 'lowering']);
        const facts = [
            { coordination: { centre: 'Cu', partner: 'F', window: [1.7, 2.52], cn: 6, whole_shell: true, confirmed: false } },
            { cluster: { element: 'Cu', members: 4, distorted: false, static: false, breathing: false, rigid: false, split: false } },
            { modulation: { elements: ['Cu'], ordered: [] } },
            { direction: { subject: 'Cu', significant: false, peak_sigma: 1, map_sigma: null } },
            { groups: [{ symbol: 'Pm-3m', ...pointGroupFacts('m-3m') }] },
            { convergence: 'improving' },
            { lowering: { parent: 'Pm-3m', lower: 'P4mm', shifted: [] } }
        ];
        const text = 'Cu is four-coordinated. The Cu4 clusters breathe. Cu forms a superstructure. Cu has a preferred direction. '
            + 'Pm-3m is polar. The run has converged. The symmetry is lowered to P4mm.';
        expect(rules(text, { facts })).toEqual(['convergence', 'polar', 'direction', 'superstructure', 'cluster', 'coordination']);
    });
});
