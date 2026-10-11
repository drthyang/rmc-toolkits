// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// Replays on the COMMITTED demo run (GaTa4Se8, 250 K) with the real engines,
// the real run context and a scripted model:
//   - the checks on the answer, and the Outcome they lower: a contradicted
//     answer becomes "not verified", a clean one keeps "achieved", and
//     gemma4:26b's correct answers of 2026-10-11 raise no check;

import { readFileSync } from 'node:fs';
import { describe, expect, it, vi } from 'vitest';
import { plotDataFromText, structureFromRmc6f } from '../browserData';
import { buildCopilotTools } from '../copilotTools';
import { runCopilotTurn } from '../llm/copilot/agentLoop';
import { contextAuditFacts } from '../llm/copilot/answerAudit';
import { buildAskUserTool } from '../llm/copilot/askUser';
import { buildRunTools } from '../llm/copilot/runTools';
import { pointGroupFacts } from '../llm/context/pointGroups';
import { buildRunContext } from '../llm/context/runContext';
import { describeSymmetry, toleranceLadder } from '../symmetryModel';
import { handlePcaMessage } from '../workers/pcaKdeWorker';

const read = (name) => readFileSync(new URL(`../../public/demo/${name}`, import.meta.url), 'utf8');
const RMC6F = read('GTS_250K.rmc6f');
const structure = structureFromRmc6f({ path: 'Demo/GTS_250K.rmc6f', text: RMC6F });
const request = (kind, params) => handlePcaMessage({ kind, ...params }, async () => RMC6F);
const partialsName = 'GTS_250K_PDFpartials.csv';
const plotFiles = [{
    name: partialsName,
    plotKind: 'pdf_partials',
    plotData: plotDataFromText({ name: partialsName, plotKind: 'pdf_partials', text: read(partialsName) })
}];
const logs = ['GTS_250K-00.log', 'GTS_250K-01.log', 'GTS_250K-02.log']
    .map((name) => plotDataFromText({ name, plotKind: 'r_value', text: read(name) }));
const rValueFile = {
    name: 'χ² history',
    plotData: { chiColumn: logs[0].chiColumn, series: [{ y: logs.flatMap((log) => log.series[0].y) }] }
};
const symmetry = { ...describeSymmetry(structure, 0.2), toleranceA: 0.2, ladder: toleranceLadder(structure, 1.0) };
const context = buildRunContext({ runName: 'Demo', plotFiles, rValueFile, structure, symmetry });
const tools = [
    ...buildRunTools({ plotFiles, rValueFile, structure, symmetry }),
    ...buildCopilotTools({ request, structure, elements: Object.keys(structure.elementCounts).sort() }),
    buildAskUserTool()
];
const auditContext = contextAuditFacts(context, pointGroupFacts);

// A scripted model: each request takes the next reply ({ text, calls }), or
// a function of the request (a model that reads what it was sent).
const scriptedModel = (replies) => {
    const requests = [];
    const callModel = (body) => {
        requests.push(body.messages.map((message) => ({ ...message })));
        const next = replies.shift();
        const reply = typeof next === 'function' ? next(body) : next;
        return (async function* stream() {
            if (!reply) throw new Error('script exhausted');
            if (reply.text) yield { content: reply.text };
            for (const [index, call] of (reply.calls || []).entries()) {
                yield { toolCall: { index, id: `c${requests.length}_${index}`, name: call.name, arguments: JSON.stringify(call.args) } };
            }
            yield { finishReason: reply.calls ? 'tool_calls' : 'stop' };
        }());
    };
    return { callModel, requests };
};

// One question: the events, the requests the model got, the cards asked.
const replay = async (question, replies, { answer = null, signal } = {}) => {
    const asked = [];
    const askUser = vi.fn(async (req) => {
        asked.push(req);
        const index = req.defaultIndex ?? 0;
        return typeof answer === 'function' ? answer(req) : (answer ?? { answer: req.options[index].label, source: 'option', index });
    });
    const { callModel, requests } = scriptedModel([...replies]);
    const events = [];
    for await (const event of runCopilotTurn({
        callModel,
        messages: [{ role: 'system', content: 'S' }, { role: 'user', content: question }],
        tools,
        signal,
        toolContext: { askUser },
        auditContext
    })) events.push(event);
    const ends = events.filter((event) => event.type === 'tool_end');
    const lastToolMessage = () => JSON.parse(requests.at(-1).filter((message) => message.role === 'tool').at(-1).content);
    return { events, requests, asked, done: events.at(-1), ends, lastToolMessage };
};


describe('the Outcome after the checks', () => {
    it('lowers a contradicted answer to "not verified" and keeps what the model claimed', async () => {
        const { done } = await replay('Has the run converged?', [
            { calls: [{ name: 'assess_convergence', args: {} }] },
            { text: 'The run has converged — the χ² history has reached a plateau.\nOutcome: achieved — flat history.' }
        ]);
        expect(done.audit.map((finding) => finding.rule)).toEqual(['convergence']);
        expect(done.outcome).toMatchObject({ status: 'not verified', claimed: 'achieved' });
        expect(done.outcome.note).toMatch(/^a check below contradicts the answer \(the model said "achieved": flat history\.\)/);
    });

    it('lowers the spot check\'s symmetry answer, from the context alone', async () => {
        const { done } = await replay('What breaks between 0.2 Å and 0.02 Å?', [
            { text: 'P-42₁m is chiral, and going from F-43m to P-42₁m the centrosymmetry is lost.\nOutcome: achieved — read the ladder.' }
        ]);
        expect(done.audit.map((finding) => finding.rule)).toEqual(['chirality', 'inversion-lost']);
        expect(done.outcome.status).toBe('not verified');
    });

    it('keeps "achieved" on gemma4:26b\'s correct answers of 2026-10-11', async () => {
        const convergence = await replay('Has the run converged?', [{
            text: 'The run has not converged; the dashboard\'s watchdog status is improving. The $\\ln(\\chi^2)$ of the last fit term '
                + '(X_ray_(R)1) has decreased from 0.0953 to −8.87 over the course of the run, and the change over the last 10% of '
                + 'the run is −0.11, indicating that the fit is still actively descending.\n'
                + 'Outcome: achieved — The convergence status was determined from the watchdog_status and ln(χ²) history.'
        }]);
        const ga = await replay('Do the Ga atoms have a preferred displacement direction?', [
            { calls: [{ name: 'list_sites', args: { element: 'Ga' } }] },
            { calls: [{ name: 'analyze_displacement_directions', args: { element: 'Ga' } }] },
            {
                text: 'Based on the displacement analysis, the Ga atoms do not exhibit a significant preferred displacement direction. '
                    + 'The pooled analysis for all Ga atoms shows no significant peak direction (significance of $-0.887\\sigma$) and '
                    + 'the overall displacement map is consistent with isotropy (map significance of $0.561\\sigma$). Furthermore, the '
                    + 'orientation tensor anisotropy is not significant ($1.1\\sigma$).\n'
                    + 'Outcome: achieved — The analysis of the Ga population shows no significant preferred direction or anisotropic orientation.'
            }
        ]);
        const space = await replay('What is the space group at 0.2 Å and at 0.02 Å tolerance, and what symmetry breaks between them?', [{
            text: 'At a 0.2 Å tolerance, the space group is F-43m (No. 216). At a 0.02 Å tolerance, the space group is P-42_1m. The '
                + 'symmetry breaking between these two tolerances is the loss of the F-43m symmetry. Note that the transition from F-43m '
                + 'to P-42_1m occurs at 0.031 Å, which is the noise_floor_A for this run, meaning this change is within the sampling noise.\n'
                + 'Outcome: achieved — The space groups and the symmetry breaking were identified from the symmetry ladder.'
        }]);
        for (const { done, asked } of [convergence, ga, space]) {
            expect(done.audit).toEqual([]);
            expect(done.outcome.status).toBe('achieved');
            expect(asked).toEqual([]);
        }
        // The Ga directions ran and their facts reached the check.
        expect(ga.ends.at(-1)).toMatchObject({ name: 'analyze_displacement_directions', status: 'ok' });
    });
});
