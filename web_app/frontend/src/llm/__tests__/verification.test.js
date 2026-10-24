// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// The verification (copilot/verification/) at its interface: record the
// steps, then judge the answer. The evidence gate's claims, its send-backs
// under both methods, numbers by reference and the computed Outcome — each
// case with the inputs and the expected output it had when the gate, the
// references and the verdict were three modules. The checks on the answer
// are in answerChecks.test.js, the catalogues' binding to the prompt in
// verificationCatalogue.test.js, the loop's use of the verification in
// agentLoop.test.js, and the validation's verdict words in
// validation/__tests__/verdictPhrasings.test.js.

import { describe, expect, it, vi } from 'vitest';
import { toolResult } from '../copilot/tools';
import { createVerification, DEFAULT_MAX_SEND_BACKS, isAnalysis, shownContext } from '../copilot/verification';

// A step as the agent loop records it: the call's result (tools.js
// toolResult: it ran and answered, unless said) with the tool's name and kind.
const step = (name, { interactive = false, reading = false, view = false, ...result } = {}) => ({
    name, interactive, reading, view, ...toolResult(result)
});
// A fail check, the tool's own refusal ('refused': it ran, without an error).
const REFUSAL = [{ level: 'fail', text: 'refused' }];
// A verification with `calls` recorded (a name, or { name, args, … }).
const verifying = (config, calls = []) => {
    const verification = createVerification(config);
    calls.forEach((call) => verification.record(step(typeof call === 'string' ? call : call.name, typeof call === 'string' ? {} : call), () => ({ data: call.data ?? null })));
    return verification;
};
// The claims an answer lacks the evidence for: the judged answer's evidence findings.
const unmet = (config, text, calls = []) => verifying(config, calls).judge(text).done.audit.filter((finding) => finding.rule === 'evidence');
const ids = (config, text, calls) => unmet(config, text, calls).map((finding) => finding.id);
// The send-back an answer gets while tools can still run, or null.
const sendBack = (config, text, calls = []) => verifying(config, calls).judge(text, { canSendBack: true }).sendBack ?? null;

const TOOLS = /\b(?:find_gr_peaks|list_sites|analyze_[a-z_]+|compute_bond_angles|detect_symmetry|measure_symmetry_lowering|compare_fit|read_run_data|assess_convergence|ask_user|read_skill)\b/;
const MO_CLUSTERS = { elements: ['Mo', 'O'], shell: (centre, partner) => (centre === partner ? { cn: centre === 'Mo' ? 2 : 6 } : null) };

describe('the evidence gate', () => {
    it('sends back a conclusion whose analysis did not run, once per claim', () => {
        const superstructure = sendBack({}, 'Ti forms a superstructure at (0, 0, ½).', ['list_sites']);
        expect(superstructure).toMatchObject({ id: 'order-between-cells', missing: ['analyze_modulations'] });
        expect(superstructure.nudge).toMatch(/^\(From the app, not the user\.\) Your answer concludes order between cells, but this question has not run analyze_modulations\./);
        expect(sendBack({}, 'Ti forms a superstructure.', ['analyze_modulations'])).toBeNull();
        // Sent back once, the same claim stands the next time.
        const once = createVerification();
        expect(once.judge('Ti forms a superstructure.', { canSendBack: true }).sendBack.id).toBe('order-between-cells');
        expect(once.judge('Ti forms a superstructure.', { canSendBack: true }).sendBack).toBeUndefined();
    });

    it('under method off (arm C) checks the same claims, and its send-back names no analysis and no tool', () => {
        const answers = [
            'Ti forms a superstructure at (0, 0, ½).',
            'Nb shows random two-site disorder along [001].',
            'Site 2 has a preferred displacement direction along [0 0 1].',
            'The Mo3 clusters breathe.',
            'Cu is six-coordinated, a Jahn–Teller octahedron.',
            'Ta is coherently shifted along c, lowering the symmetry to P4mm.',
            'There is no local order beyond thermal motion.'
        ];
        for (const answer of answers) {
            const sent = sendBack({ catalog: MO_CLUSTERS, method: 'off' }, answer);
            const named = sendBack({ catalog: MO_CLUSTERS }, answer);
            // The same claim, the same missing analyses (for the record), the same notice.
            expect({ ...sent, nudge: null }).toEqual({ ...named, nudge: null });
            expect(sent.nudge).toBe(`(From the app, not the user.) Your answer concludes ${sent.what}, but this conclusion needs its evidence `
                + 'from an analysis run in this question. Answer from the analysis\'s result (quote its numbers by reference). If none can run, '
                + 'say so in the answer instead of concluding.');
            expect(sent.nudge).not.toMatch(TOOLS);
            expect(named.nudge).toMatch(TOOLS);
        }
        expect(sendBack({ method: 'off' }, 'Ti forms a superstructure.', ['analyze_modulations'])).toBeNull();
    });

    it('reads each kind of conclusion as the answer checks read claims', () => {
        expect(ids({}, 'Nb shows random two-site disorder along [001].')).toEqual(['split-site']);
        expect(ids({}, 'Site 2 has a preferred displacement direction along [0 0 1].')).toEqual(['direction']);
        expect(ids({}, 'The Mo3 clusters breathe.')).toEqual(['clusters']);
        expect(ids({}, 'Cu is six-coordinated, a Jahn–Teller octahedron.')).toEqual(['coordination']);
        expect(ids({}, 'Ta is coherently shifted along c, lowering the symmetry to P4mm.')).toEqual(['lowering']);
        // Negated, hypothetical and asked are no conclusions.
        expect(ids({}, 'There is no superstructure, and no preferred direction.')).toEqual([]);
        expect(ids({}, 'I checked whether the Mo3 clusters breathe.')).toEqual([]);
        // The tool's own words, echoed, are not "no order" overall.
        expect(ids({}, 'No order that doubles the cell was found.')).toEqual([]);
    });

    it('"no order" needs the superstructure test, a per-site direction map and every finite cluster', () => {
        const missing = (calls) => unmet({ catalog: MO_CLUSTERS }, 'There is no local order beyond thermal motion.', calls)[0]?.missing ?? [];
        expect(missing([])).toEqual([
            'analyze_modulations',
            'analyze_displacement_directions on each kind of site (site_ref)',
            'analyze_clusters for Mo (its first Mo–Mo shell is a finite cluster\'s)'
        ]);
        // A pooled map is not a per-site one.
        expect(missing(['analyze_modulations', { name: 'analyze_displacement_directions', args: { element: 'Mo' } }])).toHaveLength(2);
        expect(missing(['analyze_modulations', { name: 'analyze_displacement_directions', args: { site_ref: 2 } },
            { name: 'analyze_clusters', args: { element: 'Mo' } }])).toEqual([]);
        // A call that failed is no evidence.
        expect(missing(['analyze_modulations', { name: 'analyze_displacement_directions', args: { site_ref: 2 } },
            { name: 'analyze_clusters', args: { element: 'Mo' }, resultStatus: 'error', error: 'boom' }])).toHaveLength(1);
    });

    it('distorted polyhedra: polyhedral_distortion where the run offers it, else analyze_clusters, either backing them (the coordination-shells interview, 2026-10-24)', () => {
        const answer = 'The TaSe6 octahedra are distorted.';
        // Every tool offered, polyhedral_distortion included: asked for it first, and either analysis meets the claim.
        expect(sendBack({}, answer)).toMatchObject({ id: 'polyhedra', missing: ['polyhedral_distortion'] });
        expect(sendBack({}, answer, ['polyhedral_distortion'])).toBeNull();
        expect(sendBack({}, answer, ['analyze_clusters'])).toBeNull();
        // A run that does not offer it: analyze_clusters is asked for, and meets it.
        const without = { offered: ['find_gr_peaks', 'compute_bond_angles', 'analyze_clusters', 'analyze_modulations'] };
        expect(sendBack(without, answer)).toMatchObject({ id: 'polyhedra', missing: ['analyze_clusters'] });
        expect(sendBack(without, answer, ['analyze_clusters'])).toBeNull();
        expect(sendBack({ offered: ['polyhedral_distortion'] }, answer)).toMatchObject({ id: 'polyhedra', missing: ['polyhedral_distortion'] });
        // A polyhedral distortion in words; a cluster's distortion is the clusters' claim, analyze_clusters alone.
        expect(ids({}, 'There is a strong octahedral distortion.')).toEqual(['polyhedra']);
        expect(ids({}, 'The Mo3 triangles are distorted.')).toEqual(['clusters']);
        expect(sendBack({}, 'The Mo3 clusters are distorted.', ['polyhedral_distortion'])).toMatchObject({ id: 'clusters', missing: ['analyze_clusters'] });
    });

    it('a split or flat-topped site called disorder goes back for the superstructure test, as random two-site disorder does (the displacement-anisotropy interview, 2026-10-24)', () => {
        expect(ids({}, 'The flat-topped Ta site is disordered.')).toEqual(['split-site']);
        expect(ids({}, 'Site 3 is disordered, a split site.')).toEqual(['split-site']);
        expect(sendBack({}, 'The flat-topped Ta site is disordered.')).toMatchObject({ id: 'split-site', missing: ['analyze_modulations'] });
        expect(sendBack({}, 'The flat-topped Ta site is disordered.', ['analyze_modulations'])).toBeNull();
        // Disorder alone, of a site the answer calls neither, is the disorder check's (it knows the sites), not the gate's.
        expect(ids({}, 'Site 3 is disordered.')).toEqual([]);
        expect(ids({}, 'The flat-topped Ta site is not disordered.')).toEqual([]);
    });

    it('anyOf asks for the first analysis the turn offers: a coordination number without find_gr_peaks offered', () => {
        const offered = ['compute_bond_angles'];
        expect(sendBack({ offered }, 'Each Ta has CN 6.')).toMatchObject({ id: 'coordination', missing: ['compute_bond_angles'] });
        expect(sendBack({ offered }, 'Each Ta has CN 6.', ['compute_bond_angles'])).toBeNull();
    });

    it('never asks for a tool the turn does not offer', () => {
        const offered = ['detect_symmetry', 'analyze_modulations'];
        expect(unmet({ offered }, 'The symmetry is lowered by a coherent shift.')[0].missing).toEqual(['detect_symmetry']);
        expect(unmet({ offered }, 'The Mo3 clusters breathe.')).toEqual([]);
        expect(createVerification({ offered }).claims.map((claim) => claim.id)).toEqual(['order-between-cells', 'split-site', 'lowering', 'no-order']);
        // No tool offered: no claim, so no gate and no evidence finding.
        expect(createVerification({ offered: [] }).claims).toEqual([]);
        expect(unmet({ offered: [] }, 'Ti forms a superstructure.')).toEqual([]);
    });

    it('reads conclusions about this run only: a definition needs no analysis', () => {
        const run = { context: { structure: { elements: { Sr: 8, Ti: 8, O: 24 } } } };
        expect(ids(run, 'A superstructure is order that enlarges the cell.')).toEqual([]);
        expect(ids(run, 'An octahedral coordination has six neighbours.')).toEqual([]);
        expect(ids(run, 'Ti forms a superstructure along c.')).toEqual(['order-between-cells']);
        expect(ids(run, 'The structure has a superstructure.')).toEqual(['order-between-cells']);
        expect(ids(run, 'TiO6 octahedra: Ti is six-coordinated.')).toEqual(['coordination']);
        // A verdict line names its element; "modulation" is the app's own word.
        expect(ids(run, 'VERDICT {"features": [{"kind": "modulation", "element": "Ti", "q": [0, 0, 0.5]}]}')).toEqual(['order-between-cells']);
        // "No order" is always about the run; a word that only contains a symbol is no mention of it.
        expect(ids(run, 'There is no local order.')).toEqual(['no-order']);
        expect(ids(run, 'Tiny clusters breathe in general.')).toEqual([]);
        // Without the run's elements, every sentence counts.
        expect(ids({ context: { run: 'r' } }, 'Ta4 tetrahedra breathe.')).toEqual(['clusters']);
    });

    it('reads a lowering in its usual words: a quantified shift or "(coherent)"; a ± shift and coherent scattering are not', () => {
        const run = { context: { structure: { elements: { K: 8, Ta: 8, O: 24 } } } };
        expect(ids(run, 'Every Ta is shifted by 0.07 Å along [0 0 1].')).toEqual(['lowering']);
        expect(ids(run, 'Every Ta shifted along [0 0 1] (coherent): the average structure is P4mm.')).toEqual(['lowering']);
        expect(ids(run, 'Ta is displaced ±0.07 Å at random.')).toEqual(['split-site']);
        expect(ids(run, 'The coherent scattering length of Ta weights the Ta–O partial.')).toEqual([]);
        // A ± sampling error, as the cluster check writes it, is no split site.
        expect(ids(run, 'The Ta edges are 3.03–3.05 Å (± 0.006 Å each).')).toEqual([]);
        expect(ids(run, 'Ta is off-centred ±0.2 Å along [0 0 1].')).toEqual(['split-site']);
        expect(ids(run, 'Every Ta is shifted by 0.07 Å along [0 0 1].', ['measure_symmetry_lowering'])).toEqual([]);
    });

    it('a random thermal displacement, a random distribution or local symmetry is no split site and no lowering (review 2026-10-17)', () => {
        const run = { context: { structure: { elements: { Ga: 4, Ta: 16, Se: 32 } } } };
        for (const text of [
            'Se shows random thermal displacements about its site.',
            'The Ta atoms are randomly distributed around their mean positions.',
            'Locally the Ta4 clusters have lower symmetry than the average structure.',
            'A local symmetry lowering of the Ta4 clusters.'
        ]) {
            expect(ids(run, text), text).toEqual([]);
        }
        expect(ids(run, 'The Ta shift lowers the average symmetry to P-42_1m.')).toEqual(['lowering']);
        expect(ids(run, 'Ta is randomly displaced along [111].')).toEqual(['split-site']);
    });
});

describe('judge: the fixed order', () => {
    it('decides the send-back on the answer with its references filled', () => {
        const verification = createVerification({ context: { structure: { elements: { Ti: 1, O: 3 } } } });
        const { ref } = verification.record(step('echo'), () => ({ data: { cn: 4 } }));
        expect(ref).toBe('r1');
        // Written by reference, the coordination number is a claim once filled.
        expect(verification.judge('Ti has CN = {{r1.cn}}.', { canSendBack: true }).sendBack).toMatchObject({ id: 'coordination', missing: ['find_gr_peaks'] });
    });

    it('sends back at most maxSendBacks times a turn (2), once per claim, then the answer stands with its findings', () => {
        expect(DEFAULT_MAX_SEND_BACKS).toBe(2);
        const text = 'Ti forms a superstructure. Site 2 has a preferred displacement direction. The Ti octahedra breathe.';
        const verification = createVerification();
        expect(verification.judge(text, { canSendBack: true, rounds: 0 }).sendBack.id).toBe('order-between-cells');
        expect(verification.judge(text, { canSendBack: true, rounds: 0 }).sendBack.id).toBe('direction');
        const { done } = verification.judge(text, { canSendBack: true, rounds: 1 });
        expect(done.gates).toEqual([
            { round: 0, id: 'order-between-cells', missing: ['analyze_modulations'] },
            { round: 0, id: 'direction', missing: ['analyze_displacement_directions'] }
        ]);
        expect(done.audit.map((finding) => finding.id)).toEqual(['order-between-cells', 'direction', 'clusters']);
        expect(createVerification({ maxSendBacks: 0 }).judge(text, { canSendBack: true }).done).toBeDefined();
        // No send-back unless tools can still run.
        expect(createVerification().judge(text).done.outcome.status).toBe('not verified');
    });

    it('blank text judges the standing answer (never sent back again), or gives null when there is none', () => {
        expect(createVerification().judge('   ')).toBeNull();
        expect(createVerification().judge(undefined)).toBeNull();
        const verification = createVerification();
        expect(verification.judge('Ti forms a superstructure.\nOutcome: sure.', { canSendBack: true }).sendBack.id).toBe('order-between-cells');
        const { done } = verification.judge('', { canSendBack: true, rounds: 1, stoppedBy: 'answer' });
        expect(done).toMatchObject({ type: 'done', text: 'Ti forms a superstructure.\nOutcome: sure.', stoppedBy: 'answer', rounds: 1 });
        expect(done.audit.map((finding) => finding.rule)).toEqual(['evidence']);
        expect(done.outcome.status).toBe('not verified');
    });

    it('done is the turn\'s event: the filled answer, the raw one, the findings, the references, the gates', () => {
        const verification = createVerification();
        verification.record(step('analyze_modulations'), () => ({ data: { elements: [{ element: 'Ti', zone_boundary: [['(0, 0, ½)', 498]] }] } }));
        const { done } = verification.judge('Ti orders at (0, 0, ½), {{r1.elements.0.zone_boundary.0.1}}× the median.\nOutcome: the point is ordered.',
            { rounds: 2, stoppedBy: 'answer', cut: false });
        expect(Object.keys(done)).toEqual(['type', 'text', 'rawText', 'outcome', 'rounds', 'stoppedBy', 'cut', 'audit', 'references', 'gates']);
        expect(done).toMatchObject({
            type: 'done', text: 'Ti orders at (0, 0, ½), 498× the median.\nOutcome: the point is ordered.', rounds: 2, stoppedBy: 'answer', cut: false,
            audit: [], references: { written: 1, filled: 1, unknown: [] }, gates: []
        });
        expect(done.outcome).toMatchObject({ status: 'achieved', claimed: null, body: 'Ti orders at (0, 0, ½), 498× the median.' });
        // Nothing to fill: no rawText.
        expect(Object.keys(createVerification().judge('Plain.').done)).not.toContain('rawText');
    });

    it('a reference that names nothing is a finding, and the number stays missing', () => {
        const verification = verifying({}, [{ name: 'echo', data: { value: 4, list: [1, 2, 3] } }]);
        const { done } = verification.judge('The value is {{r1.value}}, the list ends at {{r1.list.2}}, and the shift is {{r2.shift}} Å.');
        expect(done.text).toBe('The value is 4, the list ends at 3, and the shift is {{r2.shift}} Å.');
        expect(done.references).toEqual({ written: 3, filled: 2, unknown: ['r2.shift'] });
        expect(done.audit).toEqual([{ rule: 'reference', text: 'The answer cites {{r2.shift}}, which names no value this question\'s analyses returned: those numbers are missing.' }]);
        expect(done.outcome.status).toBe('not verified');
    });

    it('a {{ctx…}} reference names the context as the model read it, trimmed', () => {
        const ladder = Array.from({ length: 40 }, (_, index) => ({ from: index * 0.01, to: (index + 1) * 0.01, group: `G${index}`, note: 'x'.repeat(120) }));
        const context = { run: 'r', symmetry: { ladder } };
        const shown = shownContext(context);
        expect(createVerification({ context }).judge('{{ctx.symmetry.ladder.1.group}}').done.text).toBe(String(shown.symmetry.ladder[1].group));
        expect(shownContext(null)).toBeNull();
    });
});

describe('record', () => {
    it('names refs r1, r2… in record order, a step with an error taking none, and the ref resolves to exactly the data shown', () => {
        const verification = createVerification();
        const show = vi.fn((ref) => ({ data: { ref }, content: `message ${ref}` }));
        expect(verification.record(step('a'), show)).toEqual({ ref: 'r1', analysis: true, data: { ref: 'r1' }, content: 'message r1' });
        expect(verification.record(step('b', { resultStatus: 'error', error: 'boom' }), show)).toEqual({ ref: null, analysis: true, data: { ref: null }, content: 'message null' });
        expect(verification.record(step('c'), show)).toMatchObject({ ref: 'r2' });
        expect(show).toHaveBeenCalledTimes(3);
        expect(verification.judge('{{r1.ref}} and {{r2.ref}}').done.text).toBe('r1 and r2');
    });

    it('says whether the step is an analysis: the one definition (isAnalysis)', () => {
        const verification = createVerification();
        const analysis = (extra) => verification.record(step('x', extra), () => ({ data: null })).analysis;
        expect(analysis({})).toBe(true);
        expect(analysis({ resultStatus: 'error', error: 'bad arguments' })).toBe(true);
        // The tool's own refusal ran: an analysis (a compute_bond_angles refusal backs its claim).
        expect(analysis({ resultStatus: 'refused', checks: REFUSAL })).toBe(true);
        // A call that never ran is none.
        expect(analysis({ resultStatus: 'unknown_tool', error: 'unknown tool "x"' })).toBe(false);
        expect(analysis({ resultStatus: 'skipped', error: 'skipped: at most 4 tool calls per step' })).toBe(false);
        expect(analysis({ interactive: true })).toBe(false);
        expect(analysis({ reading: true })).toBe(false);
        expect(analysis({ view: true })).toBe(false);
        expect(isAnalysis({ executed: true })).toBe(true);
        expect(isAnalysis({})).toBe(false);
        expect(isAnalysis(null)).toBe(false);
    });

    it('drops a failed step\'s facts; facts that are no object are ignored, unknown keys kept and read by nothing', () => {
        const none = { direction: { subject: 'Ga', significant: false, peak_sigma: -0.887, map_sigma: 0.561 } };
        const text = 'The Ga atoms have a preferred displacement direction along [111].';
        // A result the constructors would refuse (an error with facts), for the verification's own guard.
        const broken = { ...step('analyze_displacement_directions', { resultStatus: 'error', error: 'boom' }), facts: none };
        const failed = createVerification();
        failed.record(broken, () => ({ data: null }));
        expect(failed.judge(text).done.audit.map((finding) => finding.rule)).toEqual(['evidence']);
        const odd = verifying({}, [{ name: 'analyze_displacement_directions', facts: 'significant' }, { name: 'x', facts: { someday: true } }]);
        expect(odd.judge(text).done.audit).toEqual([]);
        expect(verifying({}, [{ name: 'analyze_displacement_directions', facts: none }]).judge(text).done.audit.map((finding) => finding.rule)).toEqual(['direction']);
    });
});

describe('the computed Outcome', () => {
    it('is the checks\', never the model\'s grade', () => {
        const answer = 'Answer.\nOutcome: achieved — all good.';
        expect(verifying({}, ['list_sites']).judge(answer).done.outcome).toMatchObject({
            status: 'achieved', claimed: 'achieved', body: 'Answer.', computed: true
        });
        const none = { direction: { subject: 'Ga', significant: false, peak_sigma: -0.887, map_sigma: 0.561 } };
        const contradicted = verifying({}, [{ name: 'analyze_displacement_directions', facts: none }])
            .judge('The Ga atoms have a preferred displacement direction.\nOutcome: achieved — all good.').done.outcome;
        expect(contradicted).toMatchObject({ status: 'not verified', claimed: 'achieved' });
        expect(contradicted.note).toBe('a check below contradicts the answer. The model said "achieved": all good.');
    });

    it('unverified on missing evidence or an unfilled reference; not achieved when no analysis succeeded', () => {
        expect(createVerification().judge('Ti forms a superstructure.').done.outcome.note)
            .toBe('the analysis a conclusion needs did not run (analyze_modulations).');
        expect(createVerification().judge('{{r1.x}}').done.outcome.status).toBe('not verified');
        expect(verifying({}, [{ name: 'compute_bond_angles', resultStatus: 'refused', checks: REFUSAL }, { name: 'list_sites', resultStatus: 'error', error: 'boom' }])
            .judge('Nothing.').done.outcome).toMatchObject({ status: 'not achieved' });
    });

    it('partial when cut short, when a failure was not run again successfully, or when a check warned', () => {
        const outcome = (calls, options = {}) => verifying({}, calls).judge('An answer.', options).done.outcome;
        expect(outcome(['list_sites'], { stoppedBy: 'max_rounds' }).status).toBe('partial');
        expect(outcome([{ name: 'compute_bond_angles', resultStatus: 'refused', checks: REFUSAL }, 'list_sites']).note).toBe('compute_bond_angles failed.');
        // A refused window, run again in the whole shell: resolved.
        expect(outcome([{ name: 'compute_bond_angles', resultStatus: 'refused', checks: REFUSAL }, 'compute_bond_angles']).status).toBe('achieved');
        expect(outcome([{ name: 'find_gr_peaks', checks: [{ level: 'warn', text: 'a caveat' }] }]).note).toBe('find_gr_peaks warned: its result holds with the stated caveat.');
        // A question to the user or a skill read is no analysis.
        expect(outcome([{ name: 'ask_user', interactive: true }]).note).toBe('answered from the run context, and no check contradicts it.');
        expect(outcome([{ name: 'read_skill', reading: true }]).note).toBe('answered from the run context, and no check contradicts it.');
    });

    it('a cut reply is partial, a model may lower it, views and calls that never ran are no analyses (review 2026-10-17)', () => {
        expect(verifying({}, ['list_sites']).judge('B', { cut: true }).done.outcome.note).toMatch(/^the reply was cut at the length limit\./);
        const lowered = verifying({}, ['list_sites']).judge('B\nOutcome: not achieved — the run cannot tell.').done.outcome;
        expect(lowered).toMatchObject({ status: 'not achieved', claimed: 'not achieved' });
        expect(lowered.note).toBe('1 analysis ran and passed its checks, and no check contradicts the answer; the model says the answer is not achieved. the run cannot tell.');
        // Never raised: the checks' "not verified" stands.
        expect(createVerification().judge('{{r1.x}}\nOutcome: achieved').done.outcome.status).toBe('not verified');
        const view = { name: 'highlight_cards', resultStatus: 'error', error: 'no card', view: true };
        const typo = { name: 'analyse_sites', resultStatus: 'unknown_tool', error: 'unknown tool' };
        expect(verifying({}, [view]).judge('B').done.outcome.status).toBe('achieved');
        expect(verifying({}, [typo, 'list_sites']).judge('B').done.outcome.status).toBe('achieved');
    });
});

describe('errors', () => {
    it('createVerification refuses a phrasing for an unknown claim, a non-RegExp, or a non-global asserted pattern', () => {
        expect(() => createVerification({ phrasings: { superstructure: { asserted: [/x/g] } } })).toThrow(TypeError);
        expect(() => createVerification({ phrasings: { clusters: { asserted: ['breathe'] } } })).toThrow(TypeError);
        expect(() => createVerification({ phrasings: { clusters: { asserted: [/breathe/i] } } })).toThrow(/must be global/);
        expect(() => createVerification({ phrasings: { 'no-order': { literal: ['"features": []'] } } })).toThrow(TypeError);
        expect(() => createVerification({ phrasings: { clusters: { words: [/x/g] } } })).toThrow(TypeError);
        expect(() => createVerification({ method: 'none' })).toThrow(TypeError);
        expect(() => createVerification({ offered: 'analyze_clusters' })).toThrow(TypeError);
        expect(() => createVerification({ phrasings: { clusters: { asserted: [/(^|\W)(breathes)\b/g] } } })).not.toThrow();
    });

    it('one verification serves one turn: record or judge after done throws', () => {
        const verification = createVerification();
        verification.judge('Done.');
        expect(() => verification.judge('Again.')).toThrow(/one verification serves one turn/);
        expect(() => verification.record(step('x'), () => ({ data: null }))).toThrow(/one verification serves one turn/);
        // A send-back or a blank answer with none standing is no end.
        const open = createVerification();
        expect(open.judge('')).toBeNull();
        expect(open.judge('Ti forms a superstructure.', { canSendBack: true }).sendBack).toBeDefined();
        expect(() => open.record(step('analyze_modulations'), () => ({ data: null }))).not.toThrow();
    });

    it('model text never throws', () => {
        for (const text of [null, 42, '{{', '{{r1.}}', '{{ctx}}', 'Outcome:', '"features": [', 'VERDICT {']) {
            expect(() => createVerification({ context: { run: 'r' } }).judge(text)).not.toThrow();
        }
    });
});
