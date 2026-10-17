// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// One replay per Copilot skill (.claude/skills), on the COMMITTED demo run
// (GaTa4Se8, 250 K) with the real engines and a scripted model: the skill's
// right behaviour passes, and the old wrong one the spot check and the review
// recorded is caught by the firm rule in code — a refused tool call, or a
// check on the answer itself (llm/copilot/answerAudit.js).

import { readFileSync } from 'node:fs';
import { describe, expect, it } from 'vitest';
import { plotDataFromText, structureFromRmc6f } from '../browserData';
import { buildCopilotTools } from '../copilotTools';
import { runCopilotTurn } from '../llm/copilot/agentLoop';
import { auditAnswer } from '../llm/copilot/answerAudit';
import { buildEvidenceGate } from '../llm/copilot/evidence';
import { buildRunTools } from '../llm/copilot/runTools';
import { handlePcaMessage } from '../workers/pcaKdeWorker';

const read = (name) => readFileSync(new URL(`../../public/demo/${name}`, import.meta.url), 'utf8');
const RMC6F = read('GTS_250K.rmc6f');
const structure = structureFromRmc6f({ path: 'Demo/GTS_250K.rmc6f', text: RMC6F });
const request = (kind, params) => handlePcaMessage({ kind, ...params }, async () => RMC6F);
const hostTools = buildCopilotTools({ request, structure, elements: Object.keys(structure.elementCounts).sort() });

// The demo's own χ² history, its three logs end to end, as the Dashboard reads it.
const logs = ['GTS_250K-00.log', 'GTS_250K-01.log', 'GTS_250K-02.log']
    .map((name) => plotDataFromText({ name, plotKind: 'r_value', text: read(name) }));
const rValueFile = {
    name: 'χ² history',
    plotData: { chiColumn: logs[0].chiColumn, series: [{ y: logs.flatMap((log) => log.series[0].y) }] }
};
const tools = [...buildRunTools({ rValueFile }), ...hostTools];

// A scripted model: each request takes the next reply ({ text, calls }).
const scriptedModel = (replies) => {
    const requests = [];
    const callModel = (requestBody) => {
        requests.push(requestBody.messages.map((message) => ({ ...message })));
        const reply = replies.shift();
        return (async function* stream() {
            if (reply.text) yield { content: reply.text };
            for (const [index, call] of (reply.calls || []).entries()) {
                yield { toolCall: { index, id: `c${requests.length}_${index}`, name: call.name, arguments: JSON.stringify(call.args) } };
            }
            yield { finishReason: reply.calls ? 'tool_calls' : 'stop' };
        }());
    };
    return { callModel, requests };
};

const replay = async (replies, { auditContext = null } = {}) => {
    const { callModel, requests } = scriptedModel([...replies]);
    const events = [];
    for await (const event of runCopilotTurn({
        callModel, messages: [{ role: 'system', content: 'S' }, { role: 'user', content: 'Q' }], tools, auditContext,
        // The evidence gate is on, as in the app (THY-18).
        evidence: buildEvidenceGate({ offered: tools.map((tool) => tool.name) })
    })) events.push(event);
    const steps = events.filter((event) => event.type === 'tool_end');
    const toolMessages = requests.at(-1).filter((message) => message.role === 'tool').map((message) => JSON.parse(message.content));
    return { done: events.at(-1), steps, toolMessages };
};

describe('rmc-fit-assessment', () => {
    const ask = { calls: [{ name: 'assess_convergence', args: {} }] };

    it('reports the still-improving run as not converged', async () => {
        const { done, steps, toolMessages } = await replay([
            ask,
            { text: 'Not converged: ln χ² is still improving over the recent window.\nOutcome: achieved — the watchdog says improving.' }
        ]);
        expect(toolMessages[0].data.status).toBe('improving');
        expect(steps[0].status).toBe('ok');
        expect(done.audit).toEqual([]);
    });

    it('flags the review\'s first answer, "converged — a plateau", from the tool or from the context alone', async () => {
        const wrong = { text: 'The run has converged — the χ² history has reached a plateau.\nOutcome: achieved — flat history.' };
        const { done } = await replay([ask, wrong]);
        expect(done.audit.map((finding) => finding.rule)).toEqual(['convergence']);
        const fromContext = await replay([wrong], { auditContext: { convergence: 'improving', groups: [] } });
        expect(fromContext.done.audit.map((finding) => finding.rule)).toEqual(['convergence']);
    });
});

describe('coordination-shells', () => {
    const angles = (r12_min, r12_max) => ({ calls: [{ name: 'compute_bond_angles', args: { end1: 'Se', apex: 'Ta', end2: 'Se', r12_min, r12_max } }] });

    it('reads the whole TaSe6 shell against an ideal octahedron', async () => {
        const { steps, toolMessages } = await replay([
            angles(2.33, 2.97),
            { text: 'An octahedral split: 80 % cis near 89°, 20 % trans near 159°.\nOutcome: achieved — whole shell.' }
        ]);
        expect(steps[0].status).toBe('ok');
        const { cis_below_135_deg: cis, trans_above_135_deg: trans } = toolMessages[0].data.polyhedron;
        expect([cis.share, trans.share]).toEqual([0.8, 0.2]);
        expect(trans.mean_deg).toBeGreaterThan(150);
    });

    it('refuses the spot check\'s 2.3–2.7 Å window: no angles of a cut shell reach the model', async () => {
        const { steps, toolMessages } = await replay([
            angles(2.3, 2.7),
            { text: 'The octahedra are incomplete: coordination 3.75.\nOutcome: achieved — computed.' }
        ]);
        expect(steps[0].status).toBe('failed');
        expect(toolMessages[0].checks[0].text).toMatch(/^Refused: Only 47 % of the Ta atoms/);
        expect(toolMessages[0].data).not.toHaveProperty('modes');
        expect(toolMessages[0].data).not.toHaveProperty('mean_deg');
    });
});

describe('local-symmetry', () => {
    const ladder = { calls: [{ name: 'detect_symmetry', args: { tolerance_A: 0.2 } }, { name: 'detect_symmetry', args: { tolerance_A: 0.02 } }] };

    it('says what changes from F-43m to P-42₁m, with the right point-group facts', async () => {
        const { done, toolMessages } = await replay([
            ladder,
            { text: 'F-43m at 0.2 Å, P-42₁m at 0.02 Å: the three-fold axes and the F centring are lost. Neither group has an '
                + 'inversion centre, and P-42₁m is achiral.\nOutcome: achieved — both tolerances run.' }
        ]);
        expect(toolMessages.map((message) => message.data.point_group_facts)).toEqual([
            '-43m: cubic; no inversion centre; achiral (not a Sohncke group); non-polar',
            '-42m: tetragonal; no inversion centre; achiral (not a Sohncke group); non-polar'
        ]);
        expect(done.audit).toEqual([]);
    });

    it('flags the spot check\'s answer: "P-42₁m is chiral and centrosymmetry is lost"', async () => {
        const { done } = await replay([
            ladder,
            { text: 'P-42₁m is chiral, and going from F-43m to P-42₁m the centrosymmetry is lost.\nOutcome: achieved — x.' }
        ]);
        expect(done.audit.map((finding) => finding.rule)).toEqual(['chirality', 'inversion-lost']);
    });
});

describe('displacement-anisotropy', () => {
    const ga = { calls: [{ name: 'analyze_displacement_directions', args: { element: 'Ga' } }] };

    it('reports no preferred Ga direction below 3σ', async () => {
        const { done, steps } = await replay([
            ga,
            { text: 'No preferred direction: the peak is −0.89σ and the map 0.56σ.\nOutcome: achieved — not significant.' }
        ]);
        expect(steps[0].status).toBe('ok');
        expect(done.audit).toEqual([]);
    });

    it('flags a preferred direction the analysis did not find', async () => {
        const { done } = await replay([
            ga,
            { text: 'The Ga atoms have a preferred displacement direction along [1 1 1].\nOutcome: achieved — x.' }
        ]);
        expect(done.audit.map((finding) => finding.rule)).toEqual(['direction']);
        expect(done.audit[0].text).toMatch(/Ga: peak -0\.887σ, map 0\.561σ/);
    });
});

describe('the real model\'s correct answers (gemma4:26b, 2026-10-11)', () => {
    it('raise no check', () => {
        const groups = [
            { symbol: 'F-43m', point_group: '-43m', centrosymmetric: false, sohncke: false, polar: false },
            { symbol: 'P-42_1m', point_group: '-42m', centrosymmetric: false, sohncke: false, polar: false }
        ];
        expect(auditAnswer('The run has not converged; the dashboard\'s watchdog status is improving. The ln(χ²) of the last fit '
            + 'term (X_ray_(R)1) has decreased from 0.0953 to −8.87 over the course of the run, and the change over the last 10% '
            + 'of the run is −0.11, indicating that the fit is still actively descending.', { context: { convergence: 'improving', groups } })).toEqual([]);
        expect(auditAnswer('Based on the displacement analysis, the Ga atoms do not exhibit a significant preferred displacement '
            + 'direction. The pooled analysis for all Ga atoms shows no significant peak direction (significance of −0.887σ).', {
            facts: [{ direction: { subject: 'Ga', significant: false, peak_sigma: -0.887, map_sigma: 0.561 } }]
        })).toEqual([]);
        expect(auditAnswer('At a 0.2 Å tolerance, the space group is F-43m (No. 216). At a 0.02 Å tolerance, the space group is '
            + 'P-42_1m. The symmetry breaking between these two tolerances is the loss of the F-43m symmetry. Note that the '
            + 'transition from F-43m to P-42_1m occurs at 0.031 Å, which is the noise_floor_A for this run, meaning this change '
            + 'is within the sampling noise.', { context: { convergence: 'improving', groups } })).toEqual([]);
    });
});
