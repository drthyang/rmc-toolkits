// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

import { describe, expect, it } from 'vitest';
import { auditAnswer, contextAuditFacts, normalizeSymbol } from '../copilot/answerAudit';
import { pointGroupFacts } from '../context/pointGroups';

const rules = (text, options) => auditAnswer(text, options).map((finding) => finding.rule);

describe('the convergence check', () => {
    const improving = { facts: [{ convergence: 'improving' }] };

    it('flags "converged" against an improving run, from a tool or from the context', () => {
        expect(rules('The run has converged: the χ² history has reached a plateau.', improving)).toEqual(['convergence']);
        expect(rules('Yes — it is converged.', { context: { convergence: 'diverging', groups: [] } })).toEqual(['convergence']);
        expect(auditAnswer('It looks converged.', improving)[0].text).toMatch(/still improving .* it has not converged/);
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
    const groups = [
        { symbol: 'F-43m', ...pointGroupFacts('-43m') },
        { symbol: 'P-42_1m', ...pointGroupFacts('-42m') },
        { symbol: 'P2_13', ...pointGroupFacts('23') },
        { symbol: 'Pm-3m', ...pointGroupFacts('m-3m') }
    ];
    const context = { convergence: null, groups };

    it('reads a symbol however it is written', () => {
        expect(normalizeSymbol('P-42₁m')).toBe('p-421m');
        expect(normalizeSymbol('P−4 2_1 m')).toBe('p-421m');
        expect(normalizeSymbol('$P\\bar{4}2_1m$')).toBe('p-421m');
        expect(normalizeSymbol('F-43m (No. 216)')).toBe('f-43m');
    });

    it('flags the spot check\'s error: P-42₁m called chiral, F-43m said to lose its centrosymmetry', () => {
        const text = 'At 0.02 Å the group is P-42₁m, which is chiral. Going from F-43m to P-42₁m breaks the centrosymmetry.';
        const found = auditAnswer(text, { context });
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
        const lowered = {
            facts: [{ groups: [{ symbol: 'P4mm', ...pointGroupFacts('4mm') }] }],
            context: { convergence: null, groups: [{ symbol: 'Pm-3m', ...pointGroupFacts('m-3m'), parent: true }] }
        };
        expect(rules('At 0.01 Å the local structure loses its inversion centre.', lowered)).toEqual([]);
        expect(rules('Going to P4mm breaks the inversion symmetry of the average structure.', lowered)).toEqual([]);
        // With a parent that has no inversion, the claim is still caught.
        expect(rules('At 0.01 Å the local structure loses its inversion centre.', {
            facts: lowered.facts, context: { convergence: null, groups: [{ symbol: 'F-43m', ...pointGroupFacts('-43m'), parent: true }] }
        })).toEqual(['inversion-lost']);
    });

    it('does not blame the parent for a subgroup it does not know (review 2026-10-12 #14)', () => {
        const parentOnly = { context: { convergence: null, groups: [{ symbol: 'P-42_1m', ...pointGroupFacts('-42m'), parent: true }] } };
        expect(rules('The parent group is P-42₁m; its subgroup P2₁2₁2 is chiral.', parentOnly)).toEqual([]);
        // A clause that names no group is still read with its sentence's.
        expect(rules('The group is P-42₁m; it is chiral.', parentOnly)).toEqual(['chirality']);
    });

    it('takes the groups from a tool\'s facts too', () => {
        const facts = [{ groups: [{ symbol: 'Cmm2', ...pointGroupFacts('mm2') }] }];
        expect(rules('Cmm2 is centrosymmetric.', { facts })).toEqual(['centrosymmetry']);
        expect(rules('Cmm2 is polar.', { facts })).toEqual([]);
    });
});

describe('the preferred-direction check', () => {
    const none = { facts: [{ direction: { subject: 'Ga', significant: false, peak_sigma: -0.887, map_sigma: 0.561 } }] };

    it('flags a preferred direction no analysis of the question found', () => {
        const found = auditAnswer('The Ga atoms have a preferred displacement direction along [111].', none);
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

describe('contextAuditFacts', () => {
    it('reads the watchdog status and the ladder\'s groups from the run context', () => {
        const facts = contextAuditFacts({
            convergence: { watchdog_status: 'improving' },
            symmetry: { space_group: 'F-43m (No. 216)', point_group: '-43m', ladder: [{ sg: 'P-42_1m', pg: '-42m' }, { sg: 'P1', pg: '1' }] }
        }, pointGroupFacts);
        expect(facts.convergence).toBe('improving');
        expect(facts.groups.map((group) => [group.symbol, group.sohncke, group.centrosymmetric, Boolean(group.parent)])).toEqual([
            ['F-43m', false, false, true], ['P-42_1m', false, false, true], ['P1', true, false, true]
        ]);
        // With operation counts, the ladder's top is the rung with the most.
        const counted = contextAuditFacts({
            symmetry: { space_group: 'Pm-3m', point_group: 'm-3m', ladder: [{ sg: 'P4mm', pg: '4mm', n_ops: 8 }, { sg: 'Pm-3m', pg: 'm-3m', n_ops: 48 }, { sg: 'P1', pg: '1', n_ops: 1 }] }
        }, pointGroupFacts);
        expect(counted.groups.filter((group) => group.parent).map((group) => group.symbol)).toEqual(['Pm-3m', 'Pm-3m']);
        expect(contextAuditFacts(null, pointGroupFacts)).toBeNull();
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
        expect(auditAnswer('The Mo3 triangles are distorted.', breathing)[0].text).toMatch(/^The answer says the Mo3 clusters are distorted, but analyze_clusters found no shape distortion/);
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
        expect(auditAnswer('Cu is 4-coordinated.', shell)[0].text).toMatch(/^The answer gives Cu a coordination of 4, but its whole first shell holds 6 F per Cu in 1\.7–2\.52 Å/);
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
