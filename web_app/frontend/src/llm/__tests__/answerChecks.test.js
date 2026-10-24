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

    // The user's tiers (the displacement-anisotropy interview, 2026-10-24):
    // a preferred direction needs a finding; a hint called a hint passes.
    const tiered = (tier, because, sigma, enhancement) => ({ facts: [{
        direction: { subject: 'site 25', significant: sigma >= 3, anisotropic: true, peak_sigma: sigma, map_sigma: 5.1, tier, because, enhancement, frequency: 2 }
    }] });

    it('a peak at ≥ 3σ that is only a hint (not held, under 1.2×) or below 2.5σ backs no preferred direction', () => {
        const notHeld = audit('Site 25 has a preferred displacement direction along [0 0 1].', tiered('hint', 'not_held', 4.47, 2.234));
        expect(notHeld.map((finding) => finding.rule)).toEqual(['direction']);
        expect(notHeld[0].text).toBe('The answer claims a preferred displacement direction, but no directions analysis in it found one (a finding: ≥ 3σ, '
            + '≥ 1.2× isotropic, and held at the neighbouring resolution; site 25: peak 4.47σ at 2.234× isotropic, map 5.1σ — a hint, not held at the '
            + 'neighbouring resolution). A significant map test is anisotropy, which any thermal ellipsoid shows, not a preferred direction.');
        expect(audit('Site 25 prefers [0 0 1]: a preferred direction.', tiered('hint', 'small_enhancement', 3.4, 1.12))[0].text)
            .toMatch(/a hint, under 1\.2× isotropic/);
        expect(audit('Site 25 has a preferred direction.', tiered('none', 'below_hint_sigma', 1.9, 1.4))[0].text).toMatch(/— nothing to report\)/);
    });

    it('stays quiet on a finding, and on a hint the answer calls a hint', () => {
        expect(rules('Site 25 has a preferred displacement direction along [0 0 1].', tiered('finding', 'holds', 3.15, 1.95))).toEqual([]);
        for (const text of [
            'Site 25 shows a hint of a preferred direction along [0 0 1], not a finding.',
            'Site 25 tentatively has a preferred direction along [0 0 1].'
        ]) {
            expect(rules(text, tiered('hint', 'not_held', 4.47, 2.234)), text).toEqual([]);
        }
        // Facts without a tier (an earlier result) read their peak test.
        expect(rules('Site 1 has a preferred direction.', { facts: [{ direction: { subject: 'site 1', significant: true, peak_sigma: 4.2, map_sigma: 3.1 } }] })).toEqual([]);
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

    it('flags a coordination number the whole first shell contradicts, in the whole shell\'s window too', () => {
        // Stated in the whole shell's window (to its rounding, or as a cutoff), the claim is about the whole shell.
        expect(rules('Cu is four-coordinated (square planar) in 1.7–2.52 Å.', shell)).toEqual(['coordination']);
        expect(rules('In the whole first shell, 1.7–2.5 Å, each Cu has CN 4.', shell)).toEqual(['coordination']);
        expect(rules('Within 2.52 Å, Cu is four-coordinated.', shell)).toEqual(['coordination']);
        expect(rules('Cu sits in a square-planar coordination.', shell)).toEqual(['coordination']);
        // Without its window, the window check speaks too.
        expect(rules('Each Cu has CN 4.', shell)).toEqual(['coordination', 'window']);
        // "n A per B" is B's count of A (2026-10-24): the right count either way, the wrong one flagged.
        expect(rules('In 1.7–2.52 Å there are 4 F per Cu.', shell)).toEqual(['coordination']);
        expect(rules('In 1.7–2.52 Å there are 6 F per Cu and 1 Cu per F.', shell)).toEqual([]);
        expect(rules('The coordination in 1.7–2.52 Å is 6 F per Cu.', shell)).toEqual([]);
        expect(audit('Cu is 4-coordinated.', shell)[0].text).toMatch(/^The answer gives Cu a coordination of 4, but its whole first shell holds 6 F per Cu in 1\.7–2\.52 Å/);
    });

    it('stays quiet on the right count, a partial shell, a hypothetical, and a window the user chose', () => {
        for (const text of [
            'Cu is six-coordinated by F in 1.7–2.52 Å.',
            'Cu has an octahedral coordination with 4 short and 2 long bonds.',
            'Within the first peak (1.7–2.1 Å), Cu is four-coordinated.',
            'In 1.7–2.1 Å Cu is four-coordinated.',
            'Cu would be four-coordinated if only the short bonds counted.',
            'The coordination number of Cu is not 4.'
        ]) {
            expect(rules(text, shell)).toEqual([]);
        }
        // The user's window stands: the coordination rule is quiet, and the window check asks for the window.
        const chosen = { facts: [...shell.facts, { coordination: { centre: 'Cu', partner: 'F', window: [1.7, 2.1], cn: 4, whole_shell: false, confirmed: true } }] };
        expect(rules('Cu is four-coordinated.', chosen)).toEqual(['window']);
        expect(rules('In the window you chose, 1.7–2.1 Å, Cu is four-coordinated.', chosen)).toEqual([]);
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
        // A cluster-shape distribution: its edge ladder's width, its breathing sd.
        const clusters = [{ distribution: { kind: 'cluster_shape', subject: 'the Ta4 cluster-shape distribution', configurations: 1 } }];
        for (const text of ['The Ta4 edge ladder width is 0.12 Å against 0.10 Å for independent atoms.', 'Their breathing sd is 0.05 Å.']) {
            expect({ text, rules: rules(text, { facts: clusters }) }).toEqual({ text, rules: ['single-configuration'] });
        }
        expect(rules('The Ta4 clusters breathe: their size fluctuates 0.27 (relative) more than for independent atoms, z = 27.', { facts: clusters })).toEqual([]);
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
        expect(rules('Se–Ga–Se angles in 2.18–2.79 Å: mean 109.4°, sd 4.1°.', { facts: [{ coordination: { centre: 'Ga', partner: 'Se', window: [2.18, 2.79], cn: 4, whole_shell: true, confirmed: false } }] })).toEqual([]);
    });
});

describe('the window check: every CN and bond angle with its window in Å (the coordination-shells interview, 2026-10-24)', () => {
    // find_gr_peaks' Ta–Se and Ga–Se shells on the Demo, each way, as the tool's facts give them.
    const shells = [{ coordination: [
        { centre: 'Ta', partner: 'Se', window: [2.21, 3.23], cn: 6, whole_shell: true, confirmed: false },
        { centre: 'Se', partner: 'Ta', window: [2.21, 3.23], cn: 3, whole_shell: true, confirmed: false }
    ] }, { coordination: [
        { centre: 'Ga', partner: 'Se', window: [2.18, 2.79], cn: 3.99, whole_shell: true, confirmed: false },
        { centre: 'Se', partner: 'Ga', window: [2.18, 2.79], cn: 0.5, whole_shell: true, confirmed: false }
    ] }];

    it('flags a coordination number or a bond angle quoted without its window', () => {
        for (const text of [
            'Ta is 6-coordinated by Se.',
            'Each Ta has CN 6.',
            'There are 3.99 Se per Ga.',
            'The Se–Ta–Se angles split 80/20: cis near 89°, trans near 159°.',
            'The mean Se–Ga–Se angle is 109.4° ± 4.1°.',
            'Se–Ga–Se angles: mean 109.4°, sd 4.1°.',
            // A window stated for the Ga–Se pair carries Ga's numbers, not Ta's.
            'Ga has 4 Se neighbours in 2.18–2.79 Å. Ta has 6 Se neighbours.',
            'The cis angles average 89°.'
        ]) {
            expect({ text, rules: rules(text, { facts: shells }) }).toEqual({ text, rules: ['window'] });
        }
        expect(audit('Ta is 6-coordinated by Se.', { facts: shells })[0].text).toBe('The answer gives a coordination number ("6-coordinated") '
            + 'without its window in Å: every coordination number and bond angle carries the window it was counted in (this question\'s: '
            + 'Ta–Se 2.21–3.23 Å; Se–Ta 2.21–3.23 Å; Ga–Se 2.18–2.79 Å; Se–Ga 2.18–2.79 Å).');
        expect(audit('The cis angles average 89°.', { facts: shells })[0].text).toMatch(/^The answer gives a bond angle \(89°\) without its window in Å/);
    });

    it('stays quiet on a window in the sentence, one stated for the pair or the answer, an ideal value, a direction, and no analysis', () => {
        for (const text of [
            'Ta is 6-coordinated by Se in 2.21–3.23 Å.',
            'Each Ta is 6-coordinated by Se in [2.21, 3.23] Å.',
            'Between 2.21 and 3.23 Å each Ta has 6 Se.',
            'Within 3.23 Å each Ta has 6 Se neighbours.',
            'In the whole first shell (2.21–3.23 Å) each Ta has 6 Se. The Se–Ta–Se angles split 80/20: cis near 89°, trans near 159°.',
            'Window: the whole Ga–Se first shell, 2.18–2.79 Å.\nSe–Ga–Se angles: mean 109.4°, sd 4.1°.\nThe coordination is 3.99 Se per Ga.',
            'Window 2.18–2.79. The angles average 109.4°.',
            'An ideal octahedron has 12 cis angles at 90° and 3 trans at 180°.',
            'There are no near-180° angles.',
            'Ta is octahedrally coordinated by Se.',
            'The peak direction is [1 1 1] ± 19.4°, the cell radius at ν = 2.',
            'Ta would be 4-coordinated if only the first peak counted.'
        ]) {
            expect({ text, rules: rules(text, { facts: shells }) }).toEqual({ text, rules: [] });
        }
        expect(rules('Ta is 6-coordinated by Se.', {})).toEqual([]);
        expect(rules('The mean Se–Ga–Se angle is 109.4°.', { facts: [{ direction: { subject: 'site 1', significant: false, peak_sigma: 1, map_sigma: null } }] })).toEqual([]);
    });
});

describe('the subgroup check: a group–subgroup relation that does not hold (the local-symmetry interview, 2026-10-24)', () => {
    // The Demo's ladder: P-42₁m (8 operations of the cell) under F-43m (96).
    const context = {
        symmetry: {
            space_group: 'F-43m (No. 216)', point_group: '-43m', n_ops: 96,
            ladder: [{ sg: 'P-42_1m', pg: '-42m', holds_A: [0.02, 0.031], n_ops: 8 }, { sg: 'F-43m', pg: '-43m', holds_A: [0.031, 1], n_ops: 96 }]
        }
    };

    it('flags a relation the run\'s ladder or the point groups contradict: not a subgroup, the wrong index, t, k or maximal', () => {
        for (const [text, says] of [
            ['P-42₁m is a t-subgroup of F-43m.', /it is not translationengleiche: the translations drop from 4 to 1 per cell and the point group from -43m to -42m: a general subgroup, neither t nor k \(index 12 = 3 × 4\)\.$/],
            ['P-42₁m is a subgroup of F-43m of index 3.', /its index is 96\/8 = 12 operations of the cell, not 3:/],
            ['P-42₁m is a maximal subgroup of F-43m.', /a subgroup that is neither t nor k is never maximal/],
            ['F-43m is a subgroup of P-42₁m.', /on the ladder F-43m holds 96 operations of the cell and P-42_1m 8: P-42_1m is the subgroup, not F-43m\.$/],
            ['P-42₁m has F-43m as a subgroup.', /P-42_1m is the subgroup, not F-43m/],
            ['The symmetry lowers from F-43m to P4/mmm.', /P4\/mmm's point group, 4\/mmm, is not a subgroup of F-43m's, -43m: no subgroup of F-43m has it\.$/],
            ['F-43m → I4/mmm: the symmetry is lowered.', /I4\/mmm's point group, 4\/mmm, is not a subgroup of F-43m's/],
            ['I4/mmm is a t-subgroup of Fm-3m of index 2.', /a translationengleiche subgroup's index is the ratio of the point groups' orders, 48\/16 = 3, not 2\.$/],
            ['Fm-3m has F-43m as a k-subgroup.', /a klassengleiche subgroup keeps the point group, and F-43m's \(-43m\) is not Fm-3m's \(m-3m\)\.$/]
        ]) {
            const found = audit(text, { context });
            expect({ text, rules: found.map((finding) => finding.rule) }).toEqual({ text, rules: ['subgroup'] });
            expect(found[0].text, text).toMatch(says);
        }
    });

    it('stays quiet on a relation that holds, its right index and type, a negation, a question, and a ladder read upwards', () => {
        for (const text of [
            'Going from F-43m to P-42₁m lowers the symmetry.',
            'P-42₁m is a subgroup of F-43m of index 12, neither t nor k: the F centring and the three-fold axes are lost.',
            'I-42m is a t-subgroup of F-43m of index 3.',
            'Fm-3m → I4/mmm (t3) → P4/mmm (k2): the symmetry lowers twice.',
            'F-43m is a k-subgroup of Fm-3m? No: a t-subgroup of index 2.',
            'P4/mmm is not a subgroup of F-43m.',
            'Is P4/mmm a subgroup of F-43m?',
            'The ladder runs from P1 at 0.001 Å to F-43m at 0.031 Å.',
            'F-43m is a supergroup of P-42₁m.',
            'At 0.2 Å the group is F-43m; at 0.02 Å it is P-42₁m.'
        ]) {
            expect({ text, rules: rules(text, { context }) }).toEqual({ text, rules: [] });
        }
        // detect_symmetry's facts give a group its operations too.
        const detected = { facts: [{ groups: [{ symbol: 'P-42_1m', ...pointGroupFacts('-42m'), n_ops: 8, tolerance_A: 0.025 }] }], context };
        expect(rules('P-42₁m is a subgroup of F-43m of index 12.', detected)).toEqual([]);
        expect(rules('P-42₁m is a subgroup of F-43m of index 4.', detected)).toEqual(['subgroup']);
    });
});

describe('the box check: a cell doubling the box cannot hold (the local-symmetry interview, 2026-10-24)', () => {
    const box = (supercell) => ({ structure: { supercell, elements: { Ga: 4, Ta: 16, Se: 32 } } });

    it('flags a doubling along an axis the box is odd along, or any doubling in a box odd along every axis', () => {
        for (const [text, says] of [
            ['The Ta displacements double the cell along c.', /along c, but the box \(5 × 5 × 5 unit cells\) is odd along c: a doubled cell there does not fit in it/],
            ['This configuration shows a doubled cell.', /^The answer claims a cell doubling, but the box \(5 × 5 × 5 unit cells\) is odd along every axis/],
            ['In this structure the c axis is doubled.', /along c, but the box/],
            ['The Ta order is a cell doubling with q = (½, 0, 0).', /along a, but the box/]
        ]) {
            const found = audit(text, { context: box([5, 5, 5]) });
            expect({ text, rules: found.map((finding) => finding.rule) }).toEqual({ text, rules: ['cell-doubling'] });
            expect(found[0].text, text).toMatch(says);
        }
        expect(rules('The Ta displacements double the cell along a and c.', { context: box([4, 4, 5]) })).toEqual(['cell-doubling']);
    });

    it('stays quiet on a box that holds it, an answer that says the box cannot, a hypothetical, a denial, and a sentence not about the run', () => {
        for (const [text, supercell] of [
            ['The Ta displacements double the cell along c.', [10, 10, 10]],
            ['The Ta displacements double the cell along a.', [4, 5, 5]],
            ['The published low-temperature phase of GaTa4Se8 doubles the cell along c; this 5 × 5 × 5 box cannot hold it, so it can neither confirm nor rule out the doubling.', [5, 5, 5]],
            ['A cell doubling would need an even box.', [5, 5, 5]],
            ['No cell doubling: no zone-boundary point stands out.', [5, 5, 5]],
            ['A doubled cell is common in such phases.', [5, 5, 5]]
        ]) {
            expect({ text, rules: rules(text, { context: box(supercell) }) }).toEqual({ text, rules: [] });
        }
    });

    it('flags a doubling along axes no ordered zone-boundary point doubles, once analyze_modulations found one', () => {
        const at = (q) => ({ facts: [{ modulation: { elements: ['Ti'], ordered: [{ element: 'Ti', q }], general: [] } }], context: box([8, 8, 8]) });
        expect(rules('Ti doubles the cell along c (q = (0, 0, ½)).', at([0.5, 0, 0]))).toEqual(['cell-doubling']);
        expect(audit('Ti doubles the cell along c.', at([0.5, 0, 0]))[0].text).toBe('The answer claims a cell doubling along c, but analyze_modulations '
            + 'found ordered zone-boundary points only at (½, 0, 0): none doubles the cell along c.');
        expect(rules('Ti doubles the cell along a (q = (½, 0, 0)).', at([0.5, 0, 0]))).toEqual([]);
        expect(rules('Ti doubles the cell along c.', at([0.5, 0.5, 0.5]))).toEqual([]);
        // No ordered point at all is the superstructure check's.
        expect(rules('Ti doubles the cell along c.', { facts: [{ modulation: { elements: ['Ti'], ordered: [], general: [] } }], context: box([8, 8, 8]) }))
            .toEqual(['superstructure']);
    });
});

describe('the disorder check: a split or flat-topped site called disorder before the superstructure test (the displacement-anisotropy interview, 2026-10-24)', () => {
    // list_sites' facts: Ti site 3 split, O site 5 flat-topped, Sr site 1 neither.
    const sites = { sites: [
        { ref: 1, element: 'Sr', elements: ['Sr'], split: false, flat: false },
        { ref: 3, element: 'Ti', elements: ['Ti'], split: true, flat: false },
        { ref: 5, element: 'O', elements: ['O'], split: false, flat: true }
    ] };
    const tested = (...elements) => ({ modulation: { elements, ordered: [], general: [] } });

    it('flags disorder of a split or flat-topped site whose element no analyze_modulations tested', () => {
        for (const text of [
            'Site 3 is disordered over two positions along c.',
            'The Ti site shows positional disorder.',
            'Ti is randomly displaced along [0 0 1].',
            'Site 5 (O) is statically disordered.'
        ]) {
            expect({ text, rules: rules(text, { facts: [sites] }) }).toEqual({ text, rules: ['disorder'] });
        }
        expect(audit('Site 3 is disordered.', { facts: [sites] })[0].text).toBe('The answer calls site 3 (Ti, split) disorder, but analyze_modulations '
            + 'has not run for Ti in this question: ordered antiphase displacements (a superstructure) give a site the same split cloud as random '
            + 'disorder, and only the correlation between cells tells them apart.');
        // Tested for another element only: the O site is still untested.
        expect(rules('Site 5 is disordered.', { facts: [sites, tested('Ti')] })).toEqual(['disorder']);
    });

    it('stays quiet once its element was tested, on a ranking, a denial, thermal disorder, another site, and without list_sites', () => {
        for (const text of ['Site 3 is disordered.', 'The Ti site shows positional disorder.']) {
            expect(rules(text, { facts: [sites, tested('Ti')] }), text).toEqual([]);
            expect(rules(text, { facts: [sites, tested('Sr', 'Ti', 'O')] }), text).toEqual([]);
        }
        for (const text of [
            'Site 3 is the most disordered site: it is split, with delta_A 0.1 Å.',
            'Site 3 shows no disorder: the order at (0, 0, ½) explains the split.',
            'Ti shows thermal disorder about its site.',
            'Site 1 (Sr) is disordered.'
        ]) {
            expect({ text, rules: rules(text, { facts: [sites] }) }).toEqual({ text, rules: [] });
        }
        expect(rules('Site 3 is disordered.', {})).toEqual([]);
    });
});

describe('the order of the checks', () => {
    it('runs in FACTS order, and a lowering\'s shifts have no check yet', () => {
        expect(FACTS.map((row) => row.key)).toEqual(['convergence', 'groups', 'direction', 'sites', 'modulation', 'cluster', 'coordination', 'distribution', 'lowering']);
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
        expect(rules(text, { facts })).toEqual(['convergence', 'polar', 'direction', 'superstructure', 'cluster', 'coordination', 'window']);
    });
});
