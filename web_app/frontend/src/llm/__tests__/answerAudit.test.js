// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

import { describe, expect, it } from 'vitest';
import { auditAnswer, contextAuditFacts, lowerOutcome, normalizeSymbol } from '../copilot/answerAudit';
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

describe('lowerOutcome', () => {
    const finding = [{ rule: 'convergence', text: 'x' }];

    it('turns a contradicted "achieved" or "partial" into "not verified", keeping what the model said', () => {
        expect(lowerOutcome({ status: 'achieved', note: 'flat history.', body: 'B' }, finding)).toEqual({
            status: 'not verified',
            claimed: 'achieved',
            note: 'a check below contradicts the answer (the model said "achieved": flat history.).',
            body: 'B'
        });
        expect(lowerOutcome({ status: 'partial', note: null, body: 'B' }, [...finding, { rule: 'direction', text: 'y' }]).note)
            .toBe('2 checks below contradict the answer (the model said "partial").');
    });

    it('gives an answer without a verdict line one, and never raises "not achieved"', () => {
        expect(lowerOutcome(null, finding, 'The text.')).toEqual({
            status: 'not verified', claimed: null, note: 'a check below contradicts the answer.', body: 'The text.'
        });
        const failed = { status: 'not achieved', note: 'no data', body: 'B' };
        expect(lowerOutcome(failed, finding)).toBe(failed);
    });

    it('leaves an answer no check contradicts as the model judged it', () => {
        const fine = { status: 'achieved', note: 'n', body: 'B' };
        expect(lowerOutcome(fine, [])).toBe(fine);
        expect(lowerOutcome(null, [])).toBeNull();
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
