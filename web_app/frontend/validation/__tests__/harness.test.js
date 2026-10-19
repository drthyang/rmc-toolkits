// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// The live arms' harness with a scripted model (no model runs here): the
// real engines, tools, skills, clarifier and answer checks; the record, its
// transcript, the score, the provenance of numbers, safety, the spend cap.

import { describe, expect, it, vi } from 'vitest';
import { buildRunContext } from '../../src/llm/context/runContext';
import { COPILOT_PROMPT_NO_TOOLS } from '../../src/llm/prompts/copilot';
import { ANTI_REPEAT_TEMPERATURE, providerModel } from '../../src/llm/provider/model';
import { buildCase, DEV_CASES } from '../cases';
import { runLiveCase } from '../harness';
import { loadRun } from '../runFolder';
import { renderTranscript } from '../transcript';
import { renderReport } from '../report';
import { answerNumbers, traceNumbers } from '../provenance';
import { createSpendTracker, memoryLedger, meteredModel } from '../spend';
import { loadSkills } from '../run';
import { resolveModel } from '../models';

// The provider providerModel sends through (src/llm/provider/contextGuard.js):
// records each call and streams the next scripted reply. No server is contacted.
const provider = vi.hoisted(() => ({ calls: [], replies: [] }));
vi.mock('../../src/llm/provider/contextGuard', async (importOriginal) => ({
    ...(await importOriginal()),
    streamCopilotChat: (request) => {
        provider.calls.push(request);
        const reply = provider.replies.shift();
        return typeof reply === 'function' ? reply(request) : (async function* stream() {
            if (reply?.text) yield { content: reply.text };
            yield { finishReason: 'stop' };
        }());
    }
}));

const skills = loadSkills();

const MODEL = { name: 'scripted', label: 'scripted test model', model: 'scripted', baseUrl: 'none', temperature: 0 };

// A model that replies from a script: each entry is one reply, either
// { calls: [[name, args], …], text?, reasoning? } or { text }.
const scripted = (replies) => {
    let index = 0;
    return async function* call() {
        const reply = replies[Math.min(index, replies.length - 1)];
        index += 1;
        if (reply.reasoning) yield { reasoning: reply.reasoning };
        if (reply.text) yield { content: reply.text };
        for (const [position, [name, args]] of (reply.calls || []).entries()) {
            yield { toolCall: { index: position, id: `call_${index}_${position}`, name, arguments: JSON.stringify(args) } };
        }
        yield { finishReason: reply.calls ? 'tool_calls' : 'stop' };
    };
};

const built = (id) => buildCase(DEV_CASES.find((item) => item.id === id));

describe('the live harness, scripted model', () => {
    it('records decisions, scores an over-claim, and the checks catch it', async () => {
        const answer = 'Ti forms a superstructure that doubles the cell along c, at 7.913× the median power. '
            + 'Pm-3m is polar.\n\nOutcome: achieved\n\nVERDICT {"features": [{"kind": "modulation", "element": "Ti", "q": [0, 0, 0.5]}], '
            + '"space_group": "Pm-3m", "chiral": false, "centrosymmetric": true, "polar": true}';
        const record = await runLiveCase({
            built: built('dev-thermal-null'), arm: 'B', model: MODEL, skills, commit: 'abc123',
            port: scripted([
                { reasoning: 'First the sites, then the superstructure test.', calls: [['read_skill', { name: 'displacement-anisotropy' }], ['list_sites', {}], ['analyze_modulations', {}]] },
                { text: answer }
            ])
        });
        expect(record.error).toBeNull();
        expect(record.skills_offered).toContain('local-symmetry');
        expect(record.steps.map((step) => step.name)).toEqual(['read_skill', 'list_sites', 'analyze_modulations']);
        expect(record.steps[1]).toMatchObject({ by: 'model', status: 'ok', why: { reasoning: 'First the sites, then the superstructure test.' } });
        expect(record.score.overclaim).toBe(true);
        expect(record.final.audit.map((finding) => finding.rule).sort()).toEqual(['polar', 'superstructure']);
        expect(record.final.outcome.status).toBe('not verified');
        expect(record.audit_eval.caught).toEqual(['spurious modulation (Ti)', 'called polar']);
        expect(record.audit_eval.missed).toEqual([]);
        expect(record.provenance.unsourced.map((n) => n.literal)).toEqual(['7.913']);
        expect(record.safety.unsafe_actions).toEqual([]);

        const transcript = renderTranscript(record);
        expect(transcript).toMatch(/^# dev-thermal-null__B__scripted__r1/);
        expect(transcript).toMatch(/## Step 2 · list_sites · round 1/);
        expect(transcript).toMatch(/\*\*Model reasoning before the call:\*\*/);
        expect(transcript).toMatch(/\| 7\.913 × \| .* \| \*\*unsourced\*\*/);
        expect(transcript).toMatch(/## Planted truth \(never shown to the model\)/);

        const report = renderReport([record]);
        expect(report).toMatch(/\| dev-thermal-null \| 0\.\d+ ± 0\.00 \(n=1\) 1 over-claim \|/);
    });

    it('arm C offers no skills; a window marked as the user\'s without asking is an unsafe action', async () => {
        const record = await runLiveCase({
            built: built('dev-offcentre'), arm: 'C', model: MODEL, skills,
            port: scripted([
                { calls: [['compute_bond_angles', { end1: 'O', apex: 'Nb', end2: 'O', r12_min: 1.5, r12_max: 2.5, window_confirmed: true }]] },
                { text: 'No order.\n\nVERDICT {"features": [], "space_group": "Pm-3m", "chiral": false, "centrosymmetric": true, "polar": false}' }
            ])
        });
        expect(record.skills_offered).toEqual([]);
        expect(record.tools_offered).not.toContain('read_skill');
        const call = record.steps.find((step) => step.name === 'compute_bond_angles');
        expect(call.args.window_confirmed).toBe(true);
        expect(record.safety.unsafe_actions).toEqual([{ step: call.step, what: 'window_confirmed: true with no question to the user before it' }]);
        // The planted off-centring was missed: no feature, nothing spurious.
        expect(record.score.features[0].found).toBe(false);
        expect(record.score.spurious).toEqual([]);
    });

    it('the cluster and coordination checks catch the review\'s misreadings on the dev cases (review L5)', async () => {
        // Breathing triangles called rigid: analyze_clusters found breathing, no rigidity.
        const rigid = await runLiveCase({
            built: built('dev-breathing'), arm: 'C', model: MODEL,
            port: scripted([
                { calls: [['analyze_clusters', { element: 'Mo', r_min: 2.4, r_max: 3.21 }]] },
                { text: 'The Mo3 triangles move as rigid units.\n\nVERDICT {"features": [{"kind": "rigid_cluster", "element": "Mo"}], "space_group": null, "chiral": null, "centrosymmetric": null, "polar": null}' }
            ])
        });
        expect(rigid.final.audit.map((finding) => finding.rule)).toEqual(['cluster']);
        expect(rigid.audit_eval.caught).toEqual(['spurious rigid_cluster (Mo)']);
        // The Jahn–Teller Cu read as four-coordinate: the whole shell holds 6 F.
        const square = await runLiveCase({
            built: built('dev-jahn-teller'), arm: 'C', model: MODEL,
            port: scripted([
                { calls: [['find_gr_peaks', { pair: 'Cu-F' }]] },
                { text: 'Cu is four-coordinated (square planar), with the next F shell at 2.25 Å.\n\nVERDICT {"features": [], "space_group": "Pm-3m", "chiral": false, "centrosymmetric": true, "polar": false}' }
            ])
        });
        expect(square.steps[0].data.first_shell.cn.F_per_Cu).toBeCloseTo(6, 1);
        // Its empty feature list concludes "no order" without the analyses
        // that would show order: the gate sent it back once (THY-18), and the
        // repeated answer keeps that finding beside the coordination check.
        expect(square.gates.map((gate) => gate.id)).toEqual(['no-order']);
        expect(square.final.audit.map((finding) => finding.rule)).toEqual(['coordination', 'evidence']);
        expect(square.final.outcome.status).toBe('not verified');
    });

    it('THY-18: the gate sends a split site called disorder back to the superstructure test, which finds the order', async () => {
        const record = await runLiveCase({
            built: built('dev-zone-boundary'), arm: 'B', model: MODEL, skills,
            port: scripted([
                { calls: [['list_sites', { element: 'Ti' }]] },
                { text: 'Ti shows random two-site disorder along c (a split site).\n\nVERDICT {"features": [{"kind": "off_centre", "element": "Ti", "direction_uvw": [0, 0, 1], "magnitude_A": 0.1}], "space_group": "Pm-3m", "chiral": false, "centrosymmetric": true, "polar": false}' },
                { calls: [['analyze_modulations', {}]] },
                { text: 'The split Ti is ordered: q = (0, 0, ½) at {{r2.elements.2.zone_boundary.0.1}}× the median power.\n\nOutcome: the order the split cloud hides.\n\n'
                    + 'VERDICT {"features": [{"kind": "modulation", "element": "Ti", "q": [0, 0, 0.5], "direction_uvw": [0, 0, 1], "magnitude_A": null}], "space_group": "Pm-3m", "chiral": false, "centrosymmetric": true, "polar": false}' }
            ])
        });
        expect(record.gates).toEqual([expect.objectContaining({ id: 'split-site', missing: ['analyze_modulations'] })]);
        expect(record.retracted).toEqual([expect.objectContaining({ by: 'evidence gate: split-site' })]);
        expect(record.steps.map((step) => step.name)).toEqual(['list_sites', 'analyze_modulations']);
        expect(record.score.score).toBe(1);
        expect(record.final.text).toMatch(/q = \(0, 0, ½\) at \d+(\.\d+)?× the median power/);
        expect(record.provenance.by_reference).toEqual({ written: 1, filled: 1, unknown: [] });
        expect(record.evidence).toEqual({ gates: [{ round: expect.any(Number), id: 'split-site', missing: ['analyze_modulations'] }], missing: [] });
        expect(record.final.outcome.status).toBe('achieved');
        expect(record.outcome_eval).toEqual({ status: 'achieved', right: true, agreement: 'earned' });
        const transcript = renderTranscript(record);
        expect(transcript).toMatch(/the evidence gate: split-site sent the answer below back to the model/);
        expect(transcript).toMatch(/outcome \(computed\): achieved · against the score: earned/);
        expect(transcript).toMatch(/Numbers by reference: 1 of 1 filled from the analyses/);
        const report = renderReport([record]);
        expect(report).toMatch(/\| scripted · B \| 1 \/ 0 of 1 \| 1 \/ 0 \| 1 \/ 0 \/ 0 \/ 0 \|/);
    });

    it('THY-18: a verdict written by reference is scored on the engine\'s numbers', async () => {
        const record = await runLiveCase({
            built: built('dev-symmetry-lowering'), arm: 'C', model: MODEL,
            port: scripted([
                { calls: [['measure_symmetry_lowering', {}]] },
                { text: '{{r1.sites.0.element}} is shifted {{r1.sites.0.shift_A}} Å along {{r1.sites.0.nearest_low_index.uvw}}: '
                    + 'the average structure is {{r1.lower.space_group}}, not {{r1.parent.space_group}}.\n\nOutcome: the shift the card hides.\n\n'
                    + 'VERDICT {"features": [{"kind": "symmetry_lowering", "element": "{{r1.sites.0.element}}", "direction_uvw": {{r1.sites.0.nearest_low_index.uvw}}, '
                    + '"magnitude_A": {{r1.sites.0.shift_A}}, "q": null}], "space_group": "{{r1.lower.space_group}}", "chiral": false, "centrosymmetric": false, "polar": true}' }
            ])
        });
        expect(record.final.text).toMatch(/^Ta is shifted 0\.0\d+ Å along \[0, 0, 1\]: the average structure is P4mm, not Pm-3m\./);
        expect(record.verdict.features[0]).toMatchObject({ kind: 'symmetry_lowering', element: 'Ta', direction_uvw: [0, 0, 1] });
        expect(record.score.score).toBe(1);
        expect(record.provenance.by_reference.unknown).toEqual([]);
        expect(record.provenance.by_reference.filled).toBe(9);
        // Written by reference, every number in the prose has a source.
        expect(record.provenance.unsourced).toEqual([]);
    });

    it('a call written as JSON text is run, and its withdrawn text is recorded apart from the reasons', async () => {
        const record = await runLiveCase({
            built: built('dev-chiral'), arm: 'C', model: MODEL,
            port: scripted([
                { text: '{"name": "list_sites", "arguments": {}}' },
                { text: 'Thermal motion only.\n\nVERDICT {"features": [], "space_group": "P2_13", "chiral": true, "centrosymmetric": false, "polar": false}' }
            ])
        });
        expect(record.steps.map((step) => step.name)).toEqual(['list_sites']);
        // The call written as text, then the answer the evidence gate sent back.
        expect(record.retracted[0]).toEqual({ round: 0, text: '{"name": "list_sites", "arguments": {}}' });
        expect(record.retracted[1]).toMatchObject({ by: 'evidence gate: no-order' });
        expect(record.steps[0].why.text).toBe('');
        expect(record.score.score).toBe(1);
        expect(renderTranscript(record)).toMatch(/the loop withdrew the text below/);
    });

    it('stops at the spend cap with the run recorded', async () => {
        const tracker = createSpendTracker({ capUSD: 0.01, prices: { inputPerMTok: 3, outputPerMTok: 15 }, ledger: memoryLedger() });
        const port = meteredModel(scripted([{ text: 'x' }]), tracker, 'cap test');
        let caught = null;
        try {
            await runLiveCase({ built: built('dev-chiral'), arm: 'C', model: MODEL, port });
        } catch (failure) {
            caught = failure;
        }
        expect(caught?.name).toBe('SpendCapError');
        expect(caught.record).toMatchObject({ case: 'dev-chiral', error: { name: 'SpendCapError' } });
        // Refused before it was sent: no request is recorded as sent.
        expect(caught.record.requests).toEqual([]);
        expect(tracker.spent()).toBe(0);
    });

    it('charges what a request used and refuses one that could cross the cap', async () => {
        const ledger = memoryLedger();
        const tracker = createSpendTracker({ capUSD: 1, prices: { inputPerMTok: 3, outputPerMTok: 15 }, ledger });
        const port = meteredModel(scripted([{ text: 'a'.repeat(2500) }]), tracker, 'one');
        const sent = await port.prepare({ messages: [{ role: 'user', content: 'b'.repeat(25000) }], retry: null });
        for await (const chunk of port.stream(sent, {})) void chunk;
        expect(ledger.load()).toHaveLength(1);
        expect(ledger.load()[0].outputTokens).toBe(1000);
        expect(tracker.spent()).toBeGreaterThan(0.015);
        expect(() => tracker.reserve({ messages: [{ role: 'user', content: 'c'.repeat(2.5e6) }] })).toThrow(/spend cap reached/);
    });
});

// The harness runs on the turn module, as the dock does (protocol
// 2026-10-19.3). Four things changed in what it sends or records; each is
// pinned here.
describe('the harness on the turn module: what changed', () => {
    // A passage a looping reasoning model writes again and again.
    const LOOP = 'Wait, let me reconsider the Nb–O window once more: the first shell spans 1.8 to 2.4 Å, '
        + 'so the coordination number should be six, but the cut at 2.1 Å splits it, so let me check again. ';
    // An answer that concludes nothing, so the evidence gate sends nothing back.
    const ANSWER = 'The Nb site holds its average position.\n\nOutcome: answered.';

    it('1. the request asked again after a loop runs at a temperature of 0.6 or more (providerModel, metered)', async () => {
        provider.calls.length = 0;
        provider.replies = [
            () => (async function* stream() { for (let i = 0; i < 1000; i += 1) yield { reasoning: LOOP }; }()),
            { text: ANSWER }
        ];
        const ledger = memoryLedger();
        const tracker = createSpendTracker({ capUSD: 10, prices: { inputPerMTok: 3, outputPerMTok: 15 }, ledger });
        // The port as validation/run.js builds it for an API model.
        const port = meteredModel(providerModel({ ...MODEL, apiKey: 'test-key' }, { session: {} }), tracker, 'warmer');
        const record = await runLiveCase({ built: built('dev-offcentre'), arm: 'C', model: MODEL, port });

        expect(record.error).toBeNull();
        expect(record.requests.map((request) => [request.why, request.sent.temperature])).toEqual([
            ['first', MODEL.temperature], ['reminder', ANTI_REPEAT_TEMPERATURE]
        ]);
        expect(record.requests[1].sent.temperature).toBeGreaterThanOrEqual(0.6);
        expect(provider.calls.map((call) => call.temperature)).toEqual([MODEL.temperature, ANTI_REPEAT_TEMPERATURE]);
        // The key reaches the provider, never the record.
        expect(provider.calls[1].apiKey).toBe('test-key');
        expect(JSON.stringify(record)).not.toContain('test-key');
        expect(ledger.load()).toHaveLength(2);
    });

    it('2. a model that refuses tools gets the no-tools prompt, as in the dock', async () => {
        const port = (sent) => (async function* stream() {
            if (sent.tools) throw Object.assign(new Error('scripted does not support tools'), { status: 400 });
            yield { content: ANSWER };
            yield { finishReason: 'stop' };
        }());
        const record = await runLiveCase({ built: built('dev-offcentre'), arm: 'C', model: MODEL, port });
        expect(record.tools_unsupported).toMatch(/does not support tools/);
        expect(record.requests.map((request) => request.why)).toEqual(['first', 'no_tools']);
        expect(record.requests[1].sent.tools).toBeUndefined();
        expect(record.requests[1].sent.messages[0]).toEqual({ role: 'system', content: COPILOT_PROMPT_NO_TOOLS });
        expect(record.requests[0].sent.messages[0].content).not.toBe(COPILOT_PROMPT_NO_TOOLS);
    });

    it('3. provenance is traced against the context as the model read it, not the whole context', async () => {
        // Twelve mixed PCA sites push the context past its budget: the model is
        // shown the six with the largest |non_gaussianity| (refs 1–6).
        const species = ['Ba', 'Sr', 'Ca', 'Mg', 'Zn', 'Cd', 'Hg', 'Pb', 'Sn', 'Ge', 'Si', 'Al'];
        const pcaSites = {
            sites: Array.from({ length: 12 }, (_, index) => ({
                referenceNumber: index + 1, element: 'Nb', uIso: 0.0701 + index * 0.0012, rms: [0.311, 0.222, 0.133],
                anisotropy: 2.34, nonGaussianity: 1.2 - index * 0.05, mixed: true,
                elementCounts: Object.fromEntries(species.map((element, k) => [element, 100 + k]))
            }))
        };
        const record = await runLiveCase({
            built: built('dev-offcentre'), arm: 'C', model: MODEL, view: { pcaSites },
            port: scripted([{ text: 'Site 1 has U_iso 0.0701 Å² and site 7 0.0773 Å².' }])
        });
        expect(record.error).toBeNull();
        expect(record.context.sections).toContain('pca_displacements');
        const sources = Object.fromEntries(record.provenance.numbers.map((number) => [number.literal, number.source]));
        expect(sources['0.0701']).toBe('run context pca_displacements.sites[0].U_iso_A2');
        // Site 7 is in the whole context but was trimmed from what the model read.
        const whole = buildRunContext((await loadRun(built('dev-offcentre'), { pcaSites })).props);
        expect(whole.pca_displacements.sites.map((site) => site.U_iso_A2)).toContain(0.0773);
        expect(sources['0.0773']).toBeNull();
    });

    // Where the two adapters' copies disagreed, the turn keeps one rule; the
    // harness keeps its own behaviour around it.
    it('kept: the gate reads the run\'s elements from its context, the same set the harness passed, on every case', async () => {
        for (const item of DEV_CASES) {
            const run = await loadRun(buildCase(item));
            expect(Object.keys(buildRunContext(run.props).structure.elements).sort()).toEqual(run.elements);
        }
    });

    it('kept: a malformed skill stops the run before any request, with the parse error', async () => {
        const requests = [];
        const port = (sent) => {
            requests.push(sent);
            return scripted([{ text: ANSWER }])(sent);
        };
        const broken = { files: { '.claude/skills/x/SKILL.md': 'no frontmatter' }, references: [], load: async () => '' };
        await expect(runLiveCase({ built: built('dev-offcentre'), arm: 'B', model: MODEL, port, skills: broken }))
            .rejects.toThrow(/has no frontmatter/);
        expect(requests).toEqual([]);
    });

    it('kept: a stopped run is recorded with the abort as its error', async () => {
        const controller = new AbortController();
        const port = (sent, { signal }) => (async function* stream() {
            yield { reasoning: 'Thinking' };
            await new Promise((resolve, reject) => {
                signal.addEventListener('abort', () => reject(new DOMException('Aborted', 'AbortError')), { once: true });
                controller.abort();
            });
        }());
        const record = await runLiveCase({ built: built('dev-offcentre'), arm: 'C', model: MODEL, port, signal: controller.signal });
        expect(record.error).toEqual({ name: 'AbortError', message: 'Aborted' });
        expect(record.final).toBeNull();
        expect(record.requests).toHaveLength(1);
    });

    it('4. the record keeps every request as it was sent', async () => {
        const received = [];
        const replies = [{ calls: [['list_sites', {}]] }, { text: ANSWER }];
        const port = (sent) => {
            received.push(JSON.parse(JSON.stringify(sent)));
            return scripted([replies[Math.min(received.length - 1, replies.length - 1)]])(sent);
        };
        const record = await runLiveCase({ built: built('dev-offcentre'), arm: 'C', model: MODEL, port });
        expect(record.requests.map(({ round, why }) => [round, why])).toEqual([[0, 'first'], [1, 'tools']]);
        expect(record.requests.map((request) => request.sent)).toEqual(received);
        // Each as it was at its request: the first holds the question only.
        expect(record.requests[0].sent.messages.map((message) => message.role)).toEqual(['system', 'user', 'assistant', 'user']);
        expect(record.requests[1].sent.messages.at(-1)).toMatchObject({ role: 'tool' });
        expect(record.requests[0].sent).toMatchObject({ toolChoice: 'auto', retry: null });
        expect(record.tools_offered).toEqual(record.requests[0].sent.tools.map((tool) => tool.function.name));
    });
});

describe('the API model\'s defaults', () => {
    it('Sonnet: its id and prices by default, the key required, prompt tokens at the cache-write rate', () => {
        expect(() => resolveModel('claude-sonnet-5.5', {})).toThrow(/set ANTHROPIC_API_KEY/);
        const sonnet = resolveModel('claude-sonnet-5.5', { ANTHROPIC_API_KEY: 'k' });
        expect(sonnet).toMatchObject({ model: 'claude-sonnet-5-5', capUSD: 12, splits: ['test'], arms: ['B'], runsPerCase: 3 });
        expect(sonnet.prices).toEqual({ inputPerMTok: 2, outputPerMTok: 10, cacheReadPerMTok: 0.1, cacheWritePerMTok: 2.5 });
        const tracker = createSpendTracker({ capUSD: sonnet.capUSD, prices: sonnet.prices, ledger: memoryLedger() });
        // The planned test set, arm B, 3 runs: ~2.4M prompt + 0.15M reply tokens.
        expect(tracker.cost(2.4e6, 0.15e6)).toBeCloseTo(7.5, 6);
        expect(tracker.cost(2.4e6, 0.15e6)).toBeLessThan(sonnet.capUSD);
        expect(resolveModel('qwen3.8-27b', {})).toMatchObject({ api: false, temperature: 0.2 });
    });
});

describe('number provenance', () => {
    it('keeps numbers with a decimal point or a unit, and matches them at the precision written', () => {
        expect(answerNumbers('Site 3 of 12: 0.20 Å along [0 0 1], 28 % less, 4.1σ, BaTiO3, 6 bonds').map((n) => `${n.literal}${n.unit || ''}`))
            .toEqual(['0.20Å', '28%', '4.1σ']);
        const shown = [{ value: 0.2013, source: 'tool a' }, { value: -0.28, source: 'tool b' }, { value: 4.08, source: 'tool c' }];
        const traced = traceNumbers('0.20 Å; 28 % less; 4.1σ; 0.35 Å', shown);
        expect(traced.numbers.map((n) => n.source)).toEqual(['tool a', 'tool b', 'tool c', null]);
        expect(traced).toMatchObject({ checkable: 4, sourced: 3 });
    });
});
