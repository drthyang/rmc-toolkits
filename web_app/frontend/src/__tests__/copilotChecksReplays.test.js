// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// Replays on the COMMITTED demo run (GaTa4Se8, 250 K) with the real engines,
// the real run context and a scripted model:
//   - the checks on the answer, and the Outcome they lower: a contradicted
//     answer becomes "not verified", a clean one keeps "achieved", and
//     gemma4:26b's correct answers of 2026-10-11 raise no check;
//   - the questions the Workbench asks itself (llm/copilot/clarify.js), each
//     case, the clear questions that never ask, the model's own ask_user next
//     to them, and cancel.
// Each question runs through the turn module (llm/copilot/turn.js), as the
// app runs it: the turn builds the context, ask_user, the Workbench's
// questions, the evidence gate, the answer checks and the references from the
// run props. Two cases stay on the agent loop itself, each said where it is.

import { readFileSync } from 'node:fs';
import { describe, expect, it, vi } from 'vitest';
import { combineRValueFiles, plotDataFromText, structureFromRmc6f } from '../browserData';
import { buildCopilotTools } from '../copilotTools';
import { runCopilotTurn } from '../llm/copilot/agentLoop';
import { buildAskUserTool } from '../llm/copilot/askUser';
import { buildClarifier } from '../llm/copilot/clarify';
import { buildRunTools } from '../llm/copilot/runTools';
import { runTurn } from '../llm/copilot/turn';
import { createVerification } from '../llm/copilot/verification';
import { buildRunContext } from '../llm/context/runContext';
import { describeSymmetry, toleranceLadder } from '../symmetryModel';
import { inProcessEngine } from '../workers/pcaKdeWorker';

const read = (name) => readFileSync(new URL(`../../public/demo/${name}`, import.meta.url), 'utf8');
const RMC6F = read('GTS_250K.rmc6f');
const structure = structureFromRmc6f({ path: 'Demo/GTS_250K.rmc6f', text: RMC6F });
const request = inProcessEngine(RMC6F);
const partialsName = 'GTS_250K_PDFpartials.csv';
const plotFiles = [{
    name: partialsName,
    plotKind: 'pdf_partials',
    plotData: plotDataFromText({ name: partialsName, plotKind: 'pdf_partials', text: read(partialsName) })
}];
const logs = ['GTS_250K-00.log', 'GTS_250K-01.log', 'GTS_250K-02.log']
    .map((name) => ({ name, path: `Demo/${name}`, plotKind: 'r_value', plotData: plotDataFromText({ name, plotKind: 'r_value', text: read(name) }) }));
// Every χ² term and their total too (combineRValueFiles): what the watchdog,
// assess_convergence and the context judge.
const rValueFile = combineRValueFiles(logs, 'Demo/GTS_250K.rmc6f');
const symmetry = { ...describeSymmetry(structure, 0.2), toleranceA: 0.2, ladder: toleranceLadder(structure, 1.0) };
// The run props the turn reads, and the tools the app offers beside the
// turn's own ask_user: the tools over the run data and the engine tools.
const run = { runName: 'Demo', plotFiles, rValueFile, structure, symmetry };
const tools = [
    ...buildRunTools({ plotFiles, rValueFile, structure, symmetry }),
    ...buildCopilotTools({ request, structure, elements: Object.keys(structure.elementCounts).sort(), plotFiles })
];
const toolNamed = (name) => tools.find((tool) => tool.name === name);

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

// The user: answers each question with `answer` (an object or a function of
// the request), else picks the default option.
const answeringUser = (asked, answer) => vi.fn(async (req) => {
    asked.push(req);
    const index = req.defaultIndex ?? 0;
    return typeof answer === 'function' ? answer(req) : (answer ?? { answer: req.options[index].label, source: 'option', index });
});

// One question through the turn, as the app runs it: the loop's events (the
// turn's own start and request events left out), the requests the model got,
// the cards asked. The evidence gate is on (THY-18).
const TURN_ONLY = ['start', 'request'];
const replay = async (question, replies, { answer = null } = {}) => {
    const asked = [];
    const { callModel, requests } = scriptedModel([...replies]);
    const events = [];
    for await (const event of runTurn({ question, run, model: callModel, tools, askUser: answeringUser(asked, answer) })) {
        if (!TURN_ONLY.includes(event.type)) events.push(event);
    }
    const ends = events.filter((event) => event.type === 'tool_end');
    const lastToolMessage = () => JSON.parse(requests.at(-1).filter((message) => message.role === 'tool').at(-1).content);
    return { events, requests, asked, done: events.at(-1), ends, lastToolMessage };
};

// On the agent loop itself, for what the turn does not do: its ask_user comes
// with the Workbench's questions (one ask-user port serves both), and an abort
// ends a turn with a `stopped` event where the loop throws. The verification
// as the turn builds it.
const context = buildRunContext(run);
const clarifier = buildClarifier({ plotFiles, structure, symmetry }, { noiseFloor: context.symmetry?.noise_floor_A ?? null });
const loopReplay = async (question, replies, { answer = null } = {}) => {
    const asked = [];
    const loopTools = [...tools, buildAskUserTool()];
    const { callModel } = scriptedModel([...replies]);
    for await (const event of runCopilotTurn({
        model: callModel,
        messages: [{ role: 'system', content: 'S' }, { role: 'user', content: question }],
        tools: loopTools,
        toolContext: { askUser: answeringUser(asked, answer) },
        verification: createVerification({ context, catalog: clarifier.catalog, offered: loopTools.map((tool) => tool.name) })
    })) void event;
    return { asked };
};

const angles = (r12_min, r12_max, end = 'Se', apex = 'Ta') => ({
    calls: [{ name: 'compute_bond_angles', args: { end1: end, apex, end2: end, r12_min, r12_max } }]
});
const answer = (text) => ({ text: `${text}\n\nOutcome: achieved — computed.` });

describe('the Outcome after the checks', () => {
    it('lowers a contradicted answer to "not verified" and keeps what the model claimed', async () => {
        const { done } = await replay('Has the run converged?', [
            { calls: [{ name: 'assess_convergence', args: {} }] },
            { text: 'The run has converged — the χ² history has reached a plateau.\nOutcome: achieved — flat history.' }
        ]);
        expect(done.audit.map((finding) => finding.rule)).toEqual(['convergence']);
        expect(done.outcome).toMatchObject({ status: 'not verified', claimed: 'achieved' });
        expect(done.outcome.note).toBe('a check below contradicts the answer. The model said "achieved": flat history.');
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

describe('the questions the Workbench asks itself', () => {
    it('(a) an element only: which bond angles around Ta, then the model gets the choice', async () => {
        const run = await replay('Compute the bond angles around Ta.', [
            (body) => {
                const question = body.messages.at(-1).content;
                expect(question).toMatch(/^Compute the bond angles around Ta\.\n\n\(The Workbench asked: "Which bond angles around Ta\?" — the user chose Se–Ta–Se: Se around Ta in the whole first Ta–Se shell, 2\.21–3\.23 Å \(CN 6\)\. Use exactly this\.\)$/);
                return angles(2.21, 3.23);
            },
            answer('Octahedral.')
        ]);
        expect(run.asked).toHaveLength(1);
        expect(run.asked[0].options.map((option) => option.label)).toEqual(['Se–Ta–Se', 'Ta–Ta–Ta']);
        expect(run.asked[0].defaultIndex).toBe(0);
        // Asked before the first model request, as an ask_user step from the app.
        expect(run.events[0]).toMatchObject({ type: 'tool_start', name: 'ask_user', interactive: true, app: true });
        expect(run.events[1]).toMatchObject({ type: 'tool_end', name: 'ask_user', status: 'ok', data: { answer: 'Se–Ta–Se' } });
        // The whole shell holds the octahedron; no second question.
        expect(run.asked).toHaveLength(1);
        expect(run.lastToolMessage().data.polyhedron.trans_above_135_deg.share).toBe(0.2);

        const other = await replay('What is the coordination number of Ta?', [
            (body) => {
                expect(body.messages.at(-1).content).toMatch(/the user chose Ta around Ta: Ta around Ta in the whole first Ta–Ta shell, 2\.65–3\.33 Å \(CN 2\.8\)/);
                return answer('2.8 Ta per Ta.');
            }
        ], { answer: () => ({ answer: 'Ta around Ta', source: 'option', index: 1 }) });
        expect(other.asked[0].options.map((option) => option.label)).toEqual(['Se around Ta', 'Ta around Ta', 'Every neighbour type']);
    });

    it('(a) a split shell: which Ta–Se window, and the chosen one is not asked about again', async () => {
        // The first peak alone ends at the minimum between the two Ta–Se
        // peaks (shells.js, approved 2026-10-19): 2.21–2.65 Å, CN 3.2 — no
        // longer r ± FWHM (2.33–2.71 Å, CN 3.9).
        const run = await replay('What are the Se–Ta–Se bond angles?', [
            (body) => {
                expect(body.messages.at(-1).content).toMatch(/the user chose "First peak only" for the Ta–Se bonds: the window 2\.21–2\.65 Å \(CN 3\.2 Se per Ta\)/);
                return angles(2.21, 2.65);
            },
            answer('The short bonds only.')
        ], { answer: () => ({ answer: 'First peak only', source: 'option', index: 1 }) });
        expect(run.asked).toHaveLength(1);
        expect(run.asked[0]).toMatchObject({ question: 'Which Ta–Se window should the analysis use?', defaultIndex: 0 });
        expect(run.asked[0].options.map((option) => option.label)).toEqual(['Whole first shell', 'First peak only']);
        expect(run.asked[0].options[0].description).toBe('2.21–3.23 Å, both peaks (2.52 and 2.74 Å): CN 6 Se per Ta');
        expect(run.asked[0].options[1].description).toBe('2.21–2.65 Å, the shortest bonds: CN 3.2 Se per Ta');
        // The model's 2.21–2.65 Å is the user's choice: it runs as asked.
        expect(run.lastToolMessage().data.window_ab_A).toEqual([2.21, 2.65]);
    });

    it('(b) a window that cuts the shell: the user picks the whole shell, and the analysis runs in it', async () => {
        // The spot check's question: the octahedron settles the window, so
        // nothing is asked up front; the model's 2.3–2.7 Å is caught at the call.
        const run = await replay('Compute the Se–Ta–Se bond angles and compare them with an ideal octahedron.', [
            angles(2.3, 2.7),
            answer('Octahedral.')
        ]);
        expect(run.asked).toHaveLength(1);
        expect(run.asked[0].question).toBe('The window 2.3–2.7 Å cuts the first Ta–Se shell: it holds 3.8 Se per Ta of the shell\'s 6. Which window?');
        expect(run.asked[0].options.map((option) => option.label)).toEqual(['Whole first shell', 'Keep 2.3–2.7 Å']);
        const step = run.ends.find((end) => end.name === 'compute_bond_angles');
        expect(step.args).toMatchObject({ r12_min: 2.21, r12_max: 3.23 });
        expect(step.checks.map((check) => check.level)).toEqual(['ok', 'note']);
        const message = run.lastToolMessage();
        expect(message.note).toMatch(/^the Workbench asked the user about the window: the user chose it: run in 2\.21–3\.23 Å \(Whole first shell\)\.$/);
        expect(message.data.polyhedron.ideal).toMatch(/octahedron/);
    });

    it('(b) the user may keep the window as asked', async () => {
        const run = await replay('Compute the Se–Ta–Se bond angles and compare them with an ideal octahedron.', [
            angles(2.3, 2.7),
            answer('A partial shell.')
        ], { answer: () => ({ answer: 'Keep 2.3–2.7 Å', source: 'option', index: 1 }) });
        const step = run.ends.find((end) => end.name === 'compute_bond_angles');
        // Kept: the user chose it, so the tool runs it (window_confirmed).
        expect(step.args).toMatchObject({ r12_min: 2.3, r12_max: 2.7, window_confirmed: true });
        expect(step.status).toBe('warn');
        expect(step.checks[0].text).toMatch(/Kept: the user chose this window/);
    });

    it('(c) one site or all, for a multi-site element whose sites are not cubic', async () => {
        const run = await replay('Do the Ta atoms have a preferred displacement direction?', [
            (body) => {
                expect(body.messages.at(-1).content).toMatch(/the user chose one Ta site: analyse site_ref 1 /);
                return { calls: [{ name: 'analyze_displacement_directions', args: { site_ref: 1 } }] };
            },
            answer('No significant direction at site 1.')
        ]);
        expect(run.asked[0]).toMatchObject({ question: 'Displacement directions of Ta: one site, or all 16 pooled?' });
        expect(run.asked[0].options.map((option) => option.label)).toEqual(['One Ta site', 'All 16 Ta sites pooled']);
    });

    it('never asks a clear question', async () => {
        const clear = [
            'Has the run converged?',
            'Where is the first Ta–Se g(r) peak, and is it constraint-limited?',
            'Which site is most anharmonic, and which has the largest U_iso?',
            'Do the Ga atoms have a preferred displacement direction?',
            'Compute the Se–Ta–Se bond angles and compare them with an ideal octahedron.',
            'What is the space group at 0.2 Å and at 0.02 Å tolerance, and what symmetry breaks between them?',
            'What is the space group?',
            'Compute the Ta–Ta–Ta bond angles.',
            'Compute the Se–Ta–Se angles in 2.3–2.7 Å.'
        ];
        for (const question of clear) {
            const { asked } = await replay(question, [answer('Answered.')]);
            expect({ question, asked: asked.length }).toEqual({ question, asked: 0 });
        }
        // A window the question names is the user's: it runs without a question.
        const named = await replay('Compute the Se–Ta–Se angles in 2.3–2.7 Å.', [
            { calls: [{ name: 'compute_bond_angles', args: { end1: 'Se', apex: 'Ta', end2: 'Se', r12_min: 2.3, r12_max: 2.7, window_confirmed: true } }] },
            answer('A partial shell.')
        ]);
        expect(named.asked).toEqual([]);
        expect(named.ends[0].args).toMatchObject({ r12_min: 2.3, r12_max: 2.7 });
    });

    it('asks about a cut window once per question, then applies the choice (review 2026-10-12 #4)', async () => {
        // The model repeats its 2.3–2.7 Å call every round: the user chose the
        // whole shell the first time, and is not asked again.
        const run = await replay('Compute the Se–Ta–Se bond angles and compare them with an ideal octahedron.', [
            angles(2.3, 2.7),
            angles(2.3, 2.7),
            { calls: [{ name: 'compute_bond_angles', args: { end1: 'Se', apex: 'Ta', end2: 'Se', r12_min: 2.3, r12_max: 2.7, bin_width: 2 } }] },
            answer('Octahedral.')
        ]);
        expect(run.asked).toHaveLength(1);
        const steps = run.ends.filter((end) => end.name === 'compute_bond_angles');
        expect(steps).toHaveLength(3);
        steps.forEach((step) => expect(step.args).toMatchObject({ r12_min: 2.21, r12_max: 3.23 }));
        // The repeat of the resolved call is answered from the cache; the third
        // (other arguments) gets the remembered window and says so.
        expect(steps[1].cached).toBe(true);
        expect(run.lastToolMessage().note).toMatch(/^the Workbench applied the user's choice: the user already chose 2\.21–3\.23 Å /);
        // A window kept as asked stays kept, confirmed, on the next call too.
        const kept = await replay('Compute the Se–Ta–Se bond angles and compare them with an ideal octahedron.', [
            angles(2.3, 2.7),
            { calls: [{ name: 'compute_bond_angles', args: { end1: 'Se', apex: 'Ta', end2: 'Se', r12_min: 2.3, r12_max: 2.7, bin_width: 2 } }] },
            answer('A partial shell.')
        ], { answer: () => ({ answer: 'Keep 2.3–2.7 Å', source: 'option', index: 1 }) });
        expect(kept.asked).toHaveLength(1);
        const keptSteps = kept.ends.filter((end) => end.name === 'compute_bond_angles');
        expect(keptSteps).toHaveLength(2);
        keptSteps.forEach((step) => expect(step.args).toMatchObject({ r12_min: 2.3, r12_max: 2.7, window_confirmed: true }));
    });

    it('asks nothing when the question names the shell (review 2026-10-12 #8)', async () => {
        const first = await replay('What are the Se–Ta–Se bond angles for only the first Ta–Se peak?', [
            angles(2.33, 2.71),
            answer('The short bonds only.')
        ]);
        expect(first.asked).toEqual([]);
        // The partial window is the user's: it runs (confirmed), not refused.
        expect(first.ends[0].args).toMatchObject({ r12_min: 2.33, r12_max: 2.71, window_confirmed: true });
        expect(first.ends[0].status).toBe('warn');
        for (const question of ['What are the Se–Ta–Se bond angles in the whole first shell?', 'Se–Ta–Se angles over both Ta–Se peaks?',
            'Give the Se–Ta–Se bond angles of the entire first shell.']) {
            const whole = await replay(question, [angles(2.3, 2.7), answer('Octahedral.')]);
            expect({ question, asked: whole.asked.length }).toEqual({ question, asked: 0 });
            expect(whole.ends[0].args).toMatchObject({ r12_min: 2.21, r12_max: 3.23 });
        }
    });

    it('runs a typed window on the window the question was about (review 2026-10-12 #10)', async () => {
        // The r12 window is the whole Ta–Se shell, the r23 one cuts it: the
        // question is about r23, so a typed window goes there, not to r12.
        const cutSecond = {
            calls: [{ name: 'compute_bond_angles', args: { end1: 'Se', apex: 'Ta', end2: 'Se', r12_min: 2.21, r12_max: 3.23, r23_min: 2.3, r23_max: 2.7 } }]
        };
        const run = await replay('Compute the Se–Ta–Se bond angles and compare them with an ideal octahedron.', [cutSecond, answer('Done.')],
            { answer: () => ({ answer: '2.25 – 3.2 Å', source: 'other', index: null }) });
        expect(run.asked).toHaveLength(1);
        expect(run.ends.find((end) => end.name === 'compute_bond_angles').args)
            .toMatchObject({ r12_min: 2.21, r12_max: 3.23, r23_min: 2.25, r23_max: 3.2, window_confirmed: true });
    });

    it('leaves the model\'s own ask_user working, within two questions per request', async () => {
        const ask = { name: 'ask_user', args: { question: 'Which run?', options: ['250 K', '5 K'] } };
        const run = await replay('What are the Se–Ta–Se bond angles?', [
            { calls: [ask] },
            angles(2.3, 2.7),
            answer('Done.')
        ]);
        // The Workbench's question, then the model's: two, the limit. The
        // cut window is not asked about again: the user already chose the
        // whole shell, so it runs there (review 2026-10-12 #4, #8).
        expect(run.asked.map((req) => req.question)).toEqual(['Which Ta–Se window should the analysis use?', 'Which run?']);
        const step = run.ends.find((end) => end.name === 'compute_bond_angles');
        expect(step.status).toBe('ok');
        expect(step.args.r12_max).toBeGreaterThan(2.9);
        // Without the clarifier, the model's ask_user works alone (the loop:
        // the turn offers ask_user only with the Workbench's questions).
        const alone = await loopReplay('Compute the bond angles around Ta.', [{ calls: [ask] }, answer('Done.')]);
        expect(alone.asked.map((req) => req.question)).toEqual(['Which run?']);
    });

    it('stops cleanly when the turn is cancelled while the Workbench\'s question waits', async () => {
        // The loop throws the abort (the turn ends with a `stopped` event).
        const controller = new AbortController();
        const { callModel, requests } = scriptedModel([answer('never')]);
        const askUser = vi.fn((req, { signal }) => new Promise((_, reject) => {
            signal.addEventListener('abort', () => reject(new DOMException('Aborted', 'AbortError')), { once: true });
            setTimeout(() => controller.abort(), 5);
        }));
        const run = (async () => {
            for await (const event of runCopilotTurn({
                model: callModel, messages: [{ role: 'user', content: 'Compute the bond angles around Ta.' }], tools,
                signal: controller.signal, toolContext: { askUser }, clarifier
            })) void event;
        })();
        await expect(run).rejects.toMatchObject({ name: 'AbortError' });
        expect(askUser).toHaveBeenCalledTimes(1);
        expect(requests).toEqual([]);
    });

    it('declares the windows of the bond-angle and cluster tools', () => {
        const windows = toolNamed('compute_bond_angles').shellWindows({ end1: 'Se', apex: 'Ta', end2: 'Se', r12_min: 2.3, r12_max: 2.7 });
        expect(windows).toHaveLength(1);
        expect(windows[0]).toMatchObject({ centre: 'Ta', partner: 'Se', window: [2.3, 2.7] });
        expect(windows[0].apply([2.21, 3.23])).toMatchObject({ r12_min: 2.21, r12_max: 3.23 });
        expect(toolNamed('analyze_clusters').shellWindows({ element: 'Ta', r_min: 2.65, r_max: 3.33 })[0])
            .toMatchObject({ centre: 'Ta', partner: 'Ta', window: [2.65, 3.33] });
    });
});

describe('THY-18 on the demo run: evidence, numbers by reference, the computed verdict', () => {
    it('sends a cluster conclusion back for its analysis, fills its numbers from it, and grades it from the checks', async () => {
        const { events, done, requests } = await replay('Do the Ta4 clusters breathe?', [
            { text: 'Yes: the Ta4 clusters breathe.\nOutcome: achieved — obvious from the structure.' },
            { calls: [{ name: 'analyze_clusters', args: { element: 'Ta', r_min: 2.65, r_max: 3.33 } }] },
            {
                text: 'The Ta4 clusters breathe: their size fluctuates {{r1.local.breathing_sd_A.change}} (relative) more than for '
                    + 'independent atoms, z = {{r1.local.breathing_sd_A.z}}.\nOutcome: achieved — analyze_clusters shows breathing.'
            }
        ]);
        expect(events.find((event) => event.type === 'gate')).toMatchObject({ id: 'clusters', missing: ['analyze_clusters'] });
        expect(requests[1].at(-1).content).toMatch(/has not run analyze_clusters/);
        expect(done.text).toMatch(/^The Ta4 clusters breathe: their size fluctuates 0\.27 \(relative\) more than for independent atoms, z = 27\./);
        expect(done.references).toEqual({ written: 2, filled: 2, unknown: [] });
        expect(done.audit).toEqual([]);
        // The cluster step warned (bonds outside the window): partial, whatever the model said.
        expect(done.outcome).toMatchObject({ status: 'partial', claimed: 'achieved' });
        expect(done.outcome.note).toMatch(/^analyze_clusters warned: its result holds with the stated caveat\. The model said "achieved": analyze_clusters shows breathing\./);
    });

    it('a coordination number from the whole TaSe6 shell, by reference: achieved by the checks', async () => {
        const { done } = await replay('How many Se surround each Ta?', [
            { calls: [{ name: 'find_gr_peaks', args: { pair: 'Ta-Se' } }] },
            { text: 'Each Ta is {{r1.first_shell.cn.Se_per_Ta}}-coordinated by Se in {{r1.first_shell.window_A}} Å.\nOutcome: the whole first shell.' }
        ]);
        expect(done.text).toMatch(/^Each Ta is 6-coordinated by Se in \[2\.2\d?, 3\.\d+\] Å\./);
        expect(done.audit).toEqual([]);
        expect(done.outcome).toMatchObject({ status: 'achieved', claimed: null });
    });

    it('a wrong number typed over a reference is caught by the checks, not by the verdict alone', async () => {
        // In the whole shell's window, as every CN carries one (2026-10-24); without it the window check speaks too.
        const { done } = await replay('How many Se surround each Ta?', [
            { calls: [{ name: 'find_gr_peaks', args: { pair: 'Ta-Se' } }] },
            { text: 'Each Ta is four-coordinated by Se in 2.21–3.23 Å.\nOutcome: achieved — counted.' }
        ]);
        expect(done.audit.map((finding) => finding.rule)).toEqual(['coordination']);
        expect(done.outcome).toMatchObject({ status: 'not verified', claimed: 'achieved' });
        const bare = await replay('How many Se surround each Ta?', [
            { calls: [{ name: 'find_gr_peaks', args: { pair: 'Ta-Se' } }] },
            { text: 'Each Ta is four-coordinated by Se.\nOutcome: achieved — counted.' }
        ]);
        expect(bare.done.audit.map((finding) => finding.rule)).toEqual(['coordination', 'window']);
    });
});
