// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// THY-18: the evidence gate, numbers by reference and the computed verdict,
// each on its own (the agent loop's use of them is in agentLoop.test.js; the
// replays on the demo run and the validation's dev cases exercise them end to
// end).

import { describe, expect, it } from 'vitest';
import { buildEvidenceGate, CLUSTER_CN_MAX } from '../copilot/evidence';
import { createReferenceStore, fillReferences, formatReferenceValue } from '../copilot/references';
import { computeOutcome } from '../copilot/verdict';

const ran = (...names) => names.map((name) => (typeof name === 'string' ? { name, args: {} } : name));

describe('the evidence gate', () => {
    const gate = buildEvidenceGate();

    it('sends back a conclusion whose analysis did not run, once per rule', () => {
        const superstructure = gate.check('Ti forms a superstructure at (0, 0, ½).', { ran: ran('list_sites') });
        expect(superstructure).toMatchObject({ id: 'order-between-cells', missing: ['analyze_modulations'] });
        expect(superstructure.nudge).toMatch(/^\(From the app, not the user\.\) Your answer concludes order between cells, but this question has not run analyze_modulations\./);
        expect(gate.check('Ti forms a superstructure.', { ran: ran('analyze_modulations') })).toBeNull();
        expect(gate.check('Ti forms a superstructure.', { ran: [], fired: new Set(['order-between-cells']) })).toBeNull();
    });

    it('under method off (arm C) checks the same rules, and its send-back names no analysis and no tool', () => {
        const catalog = { elements: ['Mo', 'O'], shell: (centre, partner) => (centre === partner ? { cn: centre === 'Mo' ? 2 : 6 } : null) };
        const full = buildEvidenceGate({ catalog });
        const off = buildEvidenceGate({ catalog, method: 'off' });
        const tools = /\b(?:find_gr_peaks|list_sites|analyze_[a-z_]+|compute_bond_angles|detect_symmetry|measure_symmetry_lowering|compare_fit|read_run_data|assess_convergence|ask_user|read_skill)\b/;
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
            const sent = off.check(answer, { ran: [] });
            const named = full.check(answer, { ran: [] });
            // The same rule, the same missing analyses (for the record), the same notice.
            expect({ ...sent, nudge: null }).toEqual({ ...named, nudge: null });
            expect(sent.nudge).toBe(`(From the app, not the user.) Your answer concludes ${sent.what}, but this conclusion needs its evidence `
                + 'from an analysis run in this question. Answer from the analysis\'s result (quote its numbers by reference). If none can run, '
                + 'say so in the answer instead of concluding.');
            expect(sent.nudge).not.toMatch(tools);
            expect(named.nudge).toMatch(tools);
        }
        expect(off.check('Ti forms a superstructure.', { ran: ran('analyze_modulations') })).toBeNull();
    });

    it('reads each kind of conclusion, prose or verdict words, as the answer checks read claims', () => {
        const rule = (text, calls = []) => gate.unmet(text, ran(...calls)).map((item) => item.id);
        expect(rule('Nb shows random two-site disorder along [001].')).toEqual(['split-site']);
        expect(rule('VERDICT {"features": [{"kind": "off_centre", "element": "Nb"}]}')).toEqual(['split-site']);
        expect(rule('Site 2 has a preferred displacement direction along [0 0 1].')).toEqual(['direction']);
        expect(rule('The Mo3 clusters breathe.')).toEqual(['clusters']);
        expect(rule('Cu is six-coordinated, a Jahn–Teller octahedron.')).toEqual(['coordination']);
        expect(rule('Ta is coherently shifted along c, lowering the symmetry to P4mm.')).toEqual(['lowering']);
        // Negated, hypothetical and asked are no conclusions.
        expect(rule('There is no superstructure, and no preferred direction.')).toEqual([]);
        expect(rule('I checked whether the Mo3 clusters breathe.')).toEqual([]);
        // The tool's own words, echoed, are not "no order" overall.
        expect(rule('No order that doubles the cell was found.')).toEqual([]);
    });

    it('"no order" needs the superstructure test, a per-site direction map and every finite cluster', () => {
        const catalog = {
            elements: ['Mo', 'O'],
            shell: (centre, partner) => (centre === partner ? { cn: centre === 'Mo' ? 2 : 6 } : null)
        };
        const withClusters = buildEvidenceGate({ catalog });
        expect(withClusters.clusterElements).toEqual(['Mo']);
        expect(CLUSTER_CN_MAX).toBe(4.5);
        const missing = (calls) => withClusters.unmet('There is no local order beyond thermal motion.', ran(...calls))[0]?.missing ?? [];
        expect(missing([])).toEqual([
            'analyze_modulations',
            'analyze_displacement_directions on each kind of site (site_ref)',
            'analyze_clusters for Mo (its first Mo–Mo shell is a finite cluster\'s)'
        ]);
        // A pooled map is not a per-site one.
        expect(missing(['analyze_modulations', { name: 'analyze_displacement_directions', args: { element: 'Mo' } }])).toHaveLength(2);
        expect(missing(['analyze_modulations', { name: 'analyze_displacement_directions', args: { site_ref: 2 } },
            { name: 'analyze_clusters', args: { element: 'Mo' } }])).toEqual([]);
        // The verdict's empty feature list says the same.
        expect(withClusters.unmet('VERDICT {"features": [], "space_group": "Pm-3m"}', [])[0].id).toBe('no-order');
    });

    it('never asks for a tool the turn does not offer', () => {
        const noShift = buildEvidenceGate({ offered: ['detect_symmetry', 'analyze_modulations'] });
        expect(noShift.unmet('The symmetry is lowered by a coherent shift.', [])[0].missing).toEqual(['detect_symmetry']);
        expect(noShift.unmet('The Mo3 clusters breathe.', [])).toEqual([]);
    });
});

describe('numbers by reference', () => {
    const store = createReferenceStore({ context: { symmetry: { noise_floor_A: 0.026 } } });
    store.put('r1', { peak: { significance_sigma: 17.2, direction_uvw: [0, 0, -1] }, sites: [{ shift_A: 0.065, element: 'Ta' }] });
    store.put('r2', { local: { long_axis_within_15deg: { '<100>': { z: 12 } } } });

    it('fills each reference from the result it names, in prose and in JSON', () => {
        const out = fillReferences('Peak {{r1.peak.significance_sigma}}σ along {{ r1.peak.direction_uvw }}; {{r1.sites[0].element}} shifted '
            + '{{r1.sites.0.shift_A}} Å, floor {{ctx.symmetry.noise_floor_A}} Å; <100> z = {{r2.local.long_axis_within_15deg.<100>.z}}. '
            + 'VERDICT {"magnitude_A": {{r1.sites.0.shift_A}}, "direction_uvw": {{r1.peak.direction_uvw}}}', store);
        expect(out.text).toBe('Peak 17.2σ along [0, 0, -1]; Ta shifted 0.065 Å, floor 0.026 Å; <100> z = 12. '
            + 'VERDICT {"magnitude_A": 0.065, "direction_uvw": [0, 0, -1]}');
        expect(out).toMatchObject({ written: 8, filled: 8, unknown: [] });
        expect(JSON.parse(out.text.slice(out.text.indexOf('{')))).toEqual({ magnitude_A: 0.065, direction_uvw: [0, 0, -1] });
    });

    it('leaves a reference that names nothing as written, and anything else in braces alone', () => {
        const out = fillReferences('{{r9.x}} and {{r1.peak}} and {{r1.nope}}; $x^{{2}}$ stays.', store);
        expect(out.text).toBe('{{r9.x}} and {{r1.peak}} and {{r1.nope}}; $x^{{2}}$ stays.');
        expect(out).toMatchObject({ written: 3, filled: 0, unknown: ['r9.x', 'r1.peak', 'r1.nope'] });
        expect(fillReferences('No braces.', store)).toEqual({ text: 'No braces.', written: 0, filled: 0, unknown: [] });
        expect(formatReferenceValue(Number.NaN)).toBeNull();
        expect(formatReferenceValue({ a: 1 })).toBeNull();
    });
});

describe('the computed verdict', () => {
    const step = (name, status = 'ok', extra = {}) => ({ name, status, error: status === 'failed' && extra.error ? extra.error : null, ...extra });

    it('is the checks\', never the model\'s grade', () => {
        const model = { status: 'achieved', note: 'all good.', body: 'Answer.' };
        expect(computeOutcome({ steps: [step('list_sites')], modelOutcome: model })).toMatchObject({
            status: 'achieved', claimed: 'achieved', body: 'Answer.', computed: true
        });
        const contradicted = computeOutcome({ steps: [step('list_sites')], findings: [{ rule: 'direction', text: 'x' }], modelOutcome: model });
        expect(contradicted).toMatchObject({ status: 'not verified', claimed: 'achieved' });
        expect(contradicted.note).toBe('a check below contradicts the answer. The model said "achieved": all good.');
    });

    it('unverified on missing evidence or an unfilled reference; not achieved when no analysis succeeded', () => {
        expect(computeOutcome({ findings: [{ rule: 'evidence', missing: ['analyze_modulations'], text: 'x' }] }).note)
            .toBe('the analysis a conclusion needs did not run (analyze_modulations).');
        expect(computeOutcome({ findings: [{ rule: 'reference', text: 'x' }] }).status).toBe('not verified');
        expect(computeOutcome({ steps: [step('compute_bond_angles', 'failed'), step('list_sites', 'failed', { error: 'boom' })] }))
            .toMatchObject({ status: 'not achieved' });
    });

    it('partial when cut short, when a failure was not run again successfully, or when a check warned', () => {
        expect(computeOutcome({ steps: [step('list_sites')], stoppedBy: 'max_rounds' }).status).toBe('partial');
        expect(computeOutcome({ steps: [step('compute_bond_angles', 'failed'), step('list_sites')] }).note).toBe('compute_bond_angles failed.');
        // A refused window, run again in the whole shell: resolved.
        expect(computeOutcome({ steps: [step('compute_bond_angles', 'failed'), step('compute_bond_angles')] }).status).toBe('achieved');
        expect(computeOutcome({ steps: [step('find_gr_peaks', 'warn')] }).note).toBe('find_gr_peaks warned: its result holds with the stated caveat.');
        // A question to the user or a skill read is no analysis.
        expect(computeOutcome({ steps: [step('ask_user', 'ok', { interactive: true })] }).note).toBe('answered from the run context, and no check contradicts it.');
    });
});

describe('the turn shows the send-back and the filled answer', () => {
    it('keeps each gate\'s notice, withdraws the sent-back text, and shows the answer with its numbers filled', async () => {
        const { applyTurnEvent, newAssistantTurn } = await import('../copilot/turnState');
        let turn = newAssistantTurn(0);
        const events = [
            { type: 'round', round: 0 },
            { type: 'text', delta: 'Ti forms a superstructure.' },
            { type: 'gate', id: 'order-between-cells', notice: 'The answer concluded order between cells without analyze_modulations; the app asked for it first.' },
            { type: 'retract' },
            { type: 'round', round: 1 },
            { type: 'text', delta: 'Ordered at {{r1.x}}×.\nOutcome: ordered.' },
            {
                type: 'done', text: 'Ordered at 498×.\nOutcome: ordered.', rawText: 'Ordered at {{r1.x}}×.\nOutcome: ordered.',
                outcome: { status: 'achieved', note: '1 analysis ran.', body: 'Ordered at 498×.', claimed: null, computed: true },
                stoppedBy: 'answer', audit: [], references: { written: 1, filled: 1, unknown: [] }, gates: []
            }
        ];
        events.forEach((event) => { turn = applyTurnEvent(turn, event, 0); });
        expect(turn.gates).toEqual([{ id: 'order-between-cells', notice: expect.stringMatching(/^The answer concluded order between cells/), text: '' }]);
        expect(turn.items.filter((item) => item.kind === 'text').map((item) => item.text)).toEqual(['', 'Ordered at 498×.']);
        expect(turn.references).toEqual({ written: 1, filled: 1, unknown: [] });
    });
});

describe('the gate reads conclusions about this run only', () => {
    it('a definition names neither an element of the run nor the run, and needs no analysis', () => {
        const gate = buildEvidenceGate({ elements: ['Sr', 'Ti', 'O'] });
        expect(gate.unmet('A superstructure is order that enlarges the cell.', [])).toEqual([]);
        expect(gate.unmet('An octahedral coordination has six neighbours.', [])).toEqual([]);
        expect(gate.unmet('Ti forms a superstructure along c.', []).map((rule) => rule.id)).toEqual(['order-between-cells']);
        expect(gate.unmet('The structure has a superstructure.', []).map((rule) => rule.id)).toEqual(['order-between-cells']);
        expect(gate.unmet('TiO6 octahedra: Ti is six-coordinated.', []).map((rule) => rule.id)).toEqual(['coordination']);
        // A verdict line names its element.
        expect(gate.unmet('VERDICT {"features": [{"kind": "modulation", "element": "Ti", "q": [0, 0, 0.5]}]}', []).map((rule) => rule.id))
            .toEqual(['order-between-cells']);
        // "No order" is always about the run; a word that only contains a symbol is no mention of it.
        expect(gate.unmet('There is no local order.', []).map((rule) => rule.id)).toEqual(['no-order']);
        expect(gate.unmet('Tiny clusters breathe in general.', [])).toEqual([]);
    });
});

describe('a lowering, read in its usual words', () => {
    it('a quantified shift or "(coherent)" is a lowering; a ± shift and coherent scattering are not', () => {
        const gate = buildEvidenceGate({ elements: ['K', 'Ta', 'O'] });
        const ids = (text) => gate.unmet(text, []).map((rule) => rule.id);
        expect(ids('Every Ta is shifted by 0.07 Å along [0 0 1].')).toEqual(['lowering']);
        expect(ids('Every Ta shifted along [0 0 1] (coherent): the average structure is P4mm.')).toEqual(['lowering']);
        expect(ids('Ta is displaced ±0.07 Å at random.')).toEqual(['split-site']);
        expect(ids('The coherent scattering length of Ta weights the Ta–O partial.')).toEqual([]);
        // A ± sampling error, as the cluster check writes it, is no split site.
        expect(ids('The Ta edges are 3.03–3.05 Å (± 0.006 Å each).')).toEqual([]);
        expect(ids('Ta is off-centred ±0.2 Å along [0 0 1].')).toEqual(['split-site']);
        expect(gate.unmet('Every Ta is shifted by 0.07 Å along [0 0 1].', [{ name: 'measure_symmetry_lowering', args: {} }])).toEqual([]);
    });
});

describe('the review of the first version (2026-10-17)', () => {
    it('a random thermal displacement, a random distribution or local symmetry is no split site and no lowering', () => {
        const gate = buildEvidenceGate({ elements: ['Ga', 'Ta', 'Se'] });
        for (const text of [
            'Se shows random thermal displacements about its site.',
            'The Ta atoms are randomly distributed around their mean positions.',
            'Locally the Ta4 clusters have lower symmetry than the average structure.',
            'A local symmetry lowering of the Ta4 clusters.'
        ]) {
            expect(gate.unmet(text, []), text).toEqual([]);
        }
        expect(gate.unmet('The Ta shift lowers the average symmetry to P-42_1m.', []).map((rule) => rule.id)).toEqual(['lowering']);
        expect(gate.unmet('Ta is randomly displaced along [111].', []).map((rule) => rule.id)).toEqual(['split-site']);
        // An empty element list (a context without one) reads every sentence.
        expect(buildEvidenceGate({ elements: [] }).unmet('Ta4 tetrahedra breathe.', []).map((rule) => rule.id)).toEqual(['clusters']);
    });

    it('the verdict: a cut reply is partial, a model may lower it, views and calls that never ran are no analyses', () => {
        const ok = { name: 'list_sites', status: 'ok', executed: true };
        expect(computeOutcome({ steps: [ok], cut: true }).note).toMatch(/^the reply was cut at the length limit\./);
        const lowered = computeOutcome({ steps: [ok], modelOutcome: { status: 'not achieved', note: 'the run cannot tell.', body: 'B' } });
        expect(lowered).toMatchObject({ status: 'not achieved', claimed: 'not achieved' });
        expect(lowered.note).toBe('1 analysis ran and passed its checks, and no check contradicts the answer; the model says the answer is not achieved. the run cannot tell.');
        // Never raised: the checks' "not verified" stands.
        expect(computeOutcome({ findings: [{ rule: 'direction', text: 'x' }], modelOutcome: { status: 'achieved', note: '', body: 'B' } }).status).toBe('not verified');
        const view = { name: 'highlight_cards', status: 'failed', error: 'no card', view: true, executed: true };
        const typo = { name: 'analyse_sites', status: 'failed', error: 'unknown tool', executed: false };
        expect(computeOutcome({ steps: [view] }).status).toBe('achieved');
        expect(computeOutcome({ steps: [typo, ok] }).status).toBe('achieved');
    });

    it('a {{ctx…}} reference names the context as the model read it, trimmed', async () => {
        const { shownContext } = await import('../prompts/copilot');
        const ladder = Array.from({ length: 40 }, (_, index) => ({ from: index * 0.01, to: (index + 1) * 0.01, group: `G${index}`, note: 'x'.repeat(120) }));
        const shown = shownContext({ run: 'r', symmetry: { ladder } });
        const store = createReferenceStore({ context: shown });
        const second = fillReferences('{{ctx.symmetry.ladder.1.group}}', store).text;
        expect(second).toBe(String(shown.symmetry.ladder[1].group));
        expect(shownContext(null)).toBeNull();
    });

    it('an ended or stopped turn shows every segment with its references filled, and keeps an answer the gate withdrew', async () => {
        const { applyTurnEvent, closeTurn, newAssistantTurn } = await import('../copilot/turnState');
        let turn = newAssistantTurn(0);
        [
            { type: 'round', round: 0 },
            { type: 'tool_start', id: 't1', name: 'echo', label: 'Echo', args: {} },
            { type: 'tool_end', id: 't1', name: 'echo', status: 'ok', args: {}, data: { value: 4 }, checks: [], error: null, ms: 1, ref: 'r1' },
            { type: 'round', round: 1 },
            { type: 'text', delta: 'The value is {{r1.value}}; checking more.' }
        ].forEach((event) => { turn = applyTurnEvent(turn, event, 0); });
        const stopped = closeTurn(turn, 'stopped');
        expect(stopped.items.find((item) => item.kind === 'text' && item.text).text).toBe('The value is 4; checking more.');
        let gated = newAssistantTurn(0);
        [
            { type: 'round', round: 0 },
            { type: 'text', delta: 'Ti forms a superstructure.' },
            { type: 'gate', id: 'order-between-cells', notice: 'n', text: 'Ti forms a superstructure.' },
            { type: 'retract' }
        ].forEach((event) => { gated = applyTurnEvent(gated, event, 0); });
        expect(closeTurn(gated, 'stopped').text).toBe('Ti forms a superstructure.\n\n(Interrupted before the answer was finished.)');
    });
});
