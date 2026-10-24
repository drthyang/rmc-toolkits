// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// One replay per Copilot skill (.claude/skills), on the COMMITTED demo run
// (GaTa4Se8, 250 K) with the real engines and a scripted model: the skill's
// right behaviour passes, and the old wrong one the spot check and the review
// recorded is caught by the firm rule in code — a refused tool call, or a
// check on the answer itself (the verification, llm/copilot/verification/).

import { readFileSync } from 'node:fs';
import { describe, expect, it } from 'vitest';
import { combineRValueFiles, plotDataFromText, structureFromRmc6f } from '../browserData';
import { buildCopilotTools } from '../copilotTools';
import { buildRunTools } from '../llm/copilot/runTools';
import { buildRunContext } from '../llm/context/runContext';
import { runTurn } from '../llm/copilot/turn';
import { toolResult } from '../llm/copilot/tools';
import { createVerification } from '../llm/copilot/verification';
import { inProcessEngine } from '../workers/pcaKdeWorker';
import { buildCase, DEV_CASES } from '../../validation/cases';
import { loadRun } from '../../validation/runFolder';

const read = (name) => readFileSync(new URL(`../../public/demo/${name}`, import.meta.url), 'utf8');
const RMC6F = read('GTS_250K.rmc6f');
const structure = structureFromRmc6f({ path: 'Demo/GTS_250K.rmc6f', text: RMC6F });
const request = inProcessEngine(RMC6F);
// The demo's partials, as the Dashboard reads them: the shell catalogue
// compute_bond_angles judges a window against (the rule the Workbench asks on).
const partialsName = 'GTS_250K_PDFpartials.csv';
const plotFiles = [{
    name: partialsName,
    plotKind: 'pdf_partials',
    plotData: plotDataFromText({ name: partialsName, plotKind: 'pdf_partials', text: read(partialsName) })
}];
const hostTools = buildCopilotTools({ request, structure, elements: Object.keys(structure.elementCounts).sort(), plotFiles });

// The demo's own χ² history, its three logs end to end, as the Dashboard reads it.
const logs = ['GTS_250K-00.log', 'GTS_250K-01.log', 'GTS_250K-02.log']
    .map((name) => ({ name, path: `Demo/${name}`, plotKind: 'r_value', plotData: plotDataFromText({ name, plotKind: 'r_value', text: read(name) }) }));
// Every χ² term and their total too (combineRValueFiles): what the watchdog,
// assess_convergence and the context judge.
const rValueFile = combineRValueFiles(logs, 'Demo/GTS_250K.rmc6f');
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

// One question through the turn module (llm/copilot/turn.js), as the app runs
// it, on the run props `run` (by default the χ² history only): the turn builds
// the context, the evidence gate (on, THY-18) and the answer checks, whose
// context facts (the watchdog's verdict) come from the run.
const replay = async (replies, { run = { rValueFile }, offered = tools } = {}) => {
    const { callModel, requests } = scriptedModel([...replies]);
    const events = [];
    for await (const event of runTurn({ question: 'Q', run, model: callModel, tools: offered })) events.push(event);
    const steps = events.filter((event) => event.type === 'tool_end');
    const toolMessages = requests.at(-1).filter((message) => message.role === 'tool').map((message) => JSON.parse(message.content));
    return { done: events.at(-1), steps, toolMessages, requests };
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
        // No analysis: the run's context (the watchdog says improving) is the check.
        const fromContext = await replay([wrong]);
        expect(fromContext.done.audit.map((finding) => finding.rule)).toEqual(['convergence']);
    });
});

// The user's interview answers of 2026-10-24 (rmc-fit-assessment), each where
// a check in code can catch the wrong answer: the right answer passes and the
// old wrong one is flagged. What no check can catch is not scripted here (a
// new firm rule needs the user's word): an Rw below the data's noise taken for
// a better fit, a detailed claim about a partial under 5 % of the signal, a
// feature that differs between independent runs (the app reads one run),
// whether a subtle effect waits for the run to finish, the 10 accepted moves
// per atom past convergence, and a small effect quoted from one configuration.
describe('rmc-fit-assessment: the user\'s interview answers (2026-10-24)', () => {
    const demoRun = { rValueFile, plotFiles, structure };
    const rulesOf = (done) => done.audit.map((finding) => finding.rule);

    it('answer 1: every χ² term and the total are judged — the Demo still improves on its two data terms and the total', async () => {
        const { toolMessages } = await replay([{ calls: [{ name: 'assess_convergence', args: {} }] }, { text: 'Not converged.\nOutcome: achieved — improving.' }]);
        const { data } = toolMessages[0];
        expect(data.status).toBe('improving');
        expect(data.window_rows).toBe(227);
        // Rows 2036–2261, all after the last restart (the third log starts at row 1830).
        expect(data.terms.map((term) => [term.term, term.trend])).toEqual([
            ['F(Q)_1', 'falling'], ['Curvature #3', 'flat'], ['Curvature #5', 'flat'], ['X_ray_(R)1', 'falling']
        ]);
        expect(data.total.trend).toBe('falling');
        expect(data.zero_terms).toEqual(['Curvature #1', 'Curvature #2', 'Curvature #4', 'Curvature #6']);
        // The curvature constraint dominates the total, which falls by just over 3 × its noise.
        expect(-data.total.change_last_10pct / data.total.noise_sd).toBeGreaterThan(3);
        expect(-data.total.change_last_10pct / data.total.noise_sd).toBeLessThan(3.2);
    });

    it('answer 2, improving but slowly: a first-shell number with the "still improving" caveat passes; "converged" does not', async () => {
        const caveat = await replay([
            { calls: [{ name: 'assess_convergence', args: {} }, { name: 'compute_bond_angles', args: { end1: 'Se', apex: 'Ga', end2: 'Se', r12_min: 2.18, r12_max: 2.79 } }] },
            { text: 'The run has not converged: assess_convergence says improving (F(Q)_1, X_ray_(R)1 and the total still fall beyond their '
                + 'noise over the last tenth). With that caveat, the first shell: 3.99 Se per Ga in 2.18–2.79 Å.\nOutcome: partial — still improving.' }
        ], { run: demoRun });
        expect(rulesOf(caveat.done)).toEqual([]);
        const final = await replay([
            { calls: [{ name: 'assess_convergence', args: {} }] },
            { text: 'The run has converged, so the first shell is final: 3.99 Se per Ga.\nOutcome: achieved — converged.' }
        ], { run: demoRun });
        expect(rulesOf(final.done)).toEqual(['convergence']);
        expect(final.done.outcome.status).toBe('not verified');
    });

    it('answer 5 and rule 3(b): a distribution from the one configuration (3.5 accepted moves per atom) says so; quoted bare, it is flagged', async () => {
        const angles = { calls: [{ name: 'compute_bond_angles', args: { end1: 'Se', apex: 'Ga', end2: 'Se', r12_min: 2.18, r12_max: 2.79 } }] };
        const said = await replay([
            angles,
            { text: 'In this one configuration — lightly sampled, 3.52 accepted moves per atom, short of the 10 the method asks for — the '
                + 'Se–Ga–Se angles in 2.18–2.79 Å have mean 109.4°, sd 4.1°: indicative only, with no spread between runs to judge them by.\n'
                + 'Outcome: partial — one configuration.' }
        ], { run: demoRun });
        expect(said.toolMessages[0].data).toMatchObject({ mean_deg: 109.4, sd_deg: 4.1 });
        expect(rulesOf(said.done)).toEqual([]);
        // The README GIF's answer, as the corpus holds it: the distribution's numbers, no word of one
        // configuration (and no window: the coordination-shells interview's check, 2026-10-24).
        const bare = await replay([
            angles,
            { text: 'Se–Ga–Se angles: 23 842 angles around 4 000 Ga atoms, mean 109.4°, sd 4.1°. Mode 109.5° with share 0.94.\n'
                + 'Outcome: achieved — regular.' }
        ], { run: demoRun });
        expect(rulesOf(bare.done)).toEqual(['window', 'single-configuration']);
        expect(bare.done.audit[1].text).toMatch(/^The answer quotes the numbers of the Se–Ga–Se bond-angle distribution without saying they come from one configuration/);
        expect(bare.done.outcome.status).toBe('not verified');
        // The sampling the answer names is the context's.
        expect(buildRunContext(demoRun).configuration_optimization.accepted_moves_per_atom).toBe(3.52);
    });

    it('answer 6 and rule 3(a): order at a box-scale q is refused as a finding — Se\'s strongest mode is a long wave', async () => {
        const modes = { calls: [{ name: 'analyze_modulations', args: {} }] };
        const right = await replay([
            modes,
            { text: 'Se\'s strongest mode, q = (0.1, 0.2, -0.1) at 33× the median, is a long wave: a correlation the finite box can\'t '
                + 'represent, so it says nothing about order at its own q. No zone-boundary point stands out, so there is no superstructure.\n'
                + 'Outcome: achieved — the long waves are not order.' }
        ], { run: demoRun });
        expect(right.toolMessages[0].data.elements.find((entry) => entry.element === 'Se').strongest[0].slice(0, 2)).toEqual(['(0.1, 0.2, -0.1)', 33]);
        expect(rulesOf(right.done)).toEqual([]);
        // The long wave as order, by its q or in words; "incommensurate order"
        // is a superstructure claim too, which no ordered mode backs.
        for (const [text, rules] of [
            ['Se is ordered at q = (0.1, 0.2, -0.1): the strongest mode, 33× the median.', ['long-wave']],
            ['The Se displacements show a long-wavelength modulation, a sign of incommensurate order.', ['superstructure', 'long-wave']]
        ]) {
            const wrong = await replay([modes, { text: `${text}\nOutcome: achieved — order.` }], { run: demoRun });
            expect({ text, rules: rulesOf(wrong.done) }).toEqual({ text, rules });
        }
    });
});

describe('coordination-shells', () => {
    const angles = (r12_min, r12_max) => ({ calls: [{ name: 'compute_bond_angles', args: { end1: 'Se', apex: 'Ta', end2: 'Se', r12_min, r12_max } }] });

    it('reads the whole TaSe6 shell against an ideal octahedron', async () => {
        // In the run's one configuration, as the answer says (the single-configuration check, 2026-10-24),
        // with its window (the window check, 2026-10-24).
        const { done, steps, toolMessages } = await replay([
            angles(2.33, 2.97),
            { text: 'In this one configuration, in 2.33–2.97 Å, an octahedral split: 80 % cis near 89°, 20 % trans near 159°.\n'
                + 'Outcome: achieved — whole shell.' }
        ]);
        expect(done.audit).toEqual([]);
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
        // The CN shortfall against the whole shell, by path (the one cut rule, approved by the user 2026-10-19).
        expect(toolMessages[0].checks[0].text).toMatch(/^Refused: The window cuts the first Ta–Se shell: it holds 3\.75 Se per Ta \(shell_ab\.cn_in_window\) of the shell's 6 \(shell_ab\.cn/);
        expect(toolMessages[0].data).not.toHaveProperty('modes');
        expect(toolMessages[0].data).not.toHaveProperty('mean_deg');
    });
});

// The user's interview answers of 2026-10-24 (coordination-shells), each where
// a check in code catches the wrong answer, with the run tools over the
// partials (find_gr_peaks) offered too.
describe('coordination-shells: the user\'s interview answers (2026-10-24)', () => {
    const demoRun = { rValueFile, plotFiles, structure };
    const shellTools = [...buildRunTools(demoRun), ...hostTools];
    const rulesOf = (done) => done.audit.map((finding) => finding.rule);
    const shell = (pair = 'Ta-Se') => ({ calls: [{ name: 'find_gr_peaks', args: { pair } }] });

    it('answer 8: every CN and angle carries its window in Å — with it the answer passes, without it a warning stands under it', async () => {
        const run = { run: demoRun, offered: shellTools };
        const right = await replay([shell(), {
            text: 'In the whole Ta–Se first shell, {{r1.first_shell.window_A}} Å, each Ta has {{r1.first_shell.cn.Se_per_Ta}} Se.\n'
                + 'Outcome: achieved — the whole shell.'
        }], run);
        expect(right.done.text).toMatch(/^In the whole Ta–Se first shell, \[2\.21, 3\.23\] Å, each Ta has 6 Se\./);
        expect(rulesOf(right.done)).toEqual([]);
        for (const text of ['Each Ta is 6-coordinated by Se.', 'The Se–Ta–Se angles split 80/20: cis near 89°, trans near 159°.']) {
            const wrong = await replay([shell(), { text: `${text}\nOutcome: achieved — counted.` }], run);
            expect({ text, rules: rulesOf(wrong.done) }).toEqual({ text, rules: ['window'] });
            expect(wrong.done.audit[0].text).toMatch(/without its window in Å: every coordination number and bond angle carries the window it was counted in \(this question's: Ta–Se 2\.21–3\.23 Å; Se–Ta 2\.21–3\.23 Å\)\.$/);
            expect(wrong.done.outcome.status).toBe('not verified');
        }
    });

    it('answer 1: find_gr_peaks cross-checks the window against the average structure\'s bond lengths — every Demo shell holds them', async () => {
        const pairs = ['Ga-Se', 'Ta-Se', 'Ta-Ta', 'Se-Se'];
        const { toolMessages } = await replay([
            { calls: pairs.map((pair) => ({ name: 'find_gr_peaks', args: { pair } })) },
            { text: 'Read.\nOutcome: achieved — read.' }
        ], { run: demoRun, offered: shellTools });
        // The average structure's first shell for each pair (shortest and longest bond between the site means), inside its window.
        expect(toolMessages.map((message) => [message.data.pair, message.data.first_shell.window_A, message.data.first_shell.average_structure_d_A])).toEqual([
            ['Ga-Se', [2.18, 2.79], [2.42, 2.44]],
            ['Ta-Se', [2.21, 3.23], [2.52, 2.76]],
            ['Ta-Ta', [2.65, 3.33], [3.02, 3.07]],
            ['Se-Se', [3.02, 4.42], [3.35, 3.69]]
        ]);
        expect(toolMessages.flatMap((message) => message.checks).filter((check) => /average structure/.test(check.text))).toEqual([]);
    });

    // The title is the one its decision record names (2026-10-24-coordination-shells-interview); the
    // Demo now offers polyhedral_distortion, so both cases run: offered, then not.
    it('answer 7: distorted polyhedra need a distortion index where the run offers one, else the cluster analysis (the Demo offers no polyhedral_distortion today)', async () => {
        // Offered (the Demo's tools): the gate asks for polyhedral_distortion.
        expect(shellTools.find((tool) => tool.name === 'polyhedral_distortion')).toBeDefined();
        const offered = await replay([
            { text: 'The TaSe6 octahedra are distorted.\nOutcome: achieved — obvious.' },
            { calls: [{ name: 'polyhedral_distortion', args: { centre: 'Ta', ligand: 'Se', r_min: 2.21, r_max: 3.23 } }] },
            { text: 'In this one configuration, polyhedral_distortion reads the TaSe6 octahedra.\nOutcome: partial — one configuration.' }
        ], { run: demoRun, offered: shellTools });
        expect(offered.done.gates).toEqual([{ round: 0, id: 'polyhedra', missing: ['polyhedral_distortion'] }]);
        expect(offered.requests[1].at(-1).content).toMatch(/concludes a distortion of the coordination polyhedra, but this question has not run polyhedral_distortion\./);
        // Not offered: the cluster analysis.
        const { done, requests } = await replay([
            { text: 'The TaSe6 octahedra are distorted.\nOutcome: achieved — obvious.' },
            { calls: [{ name: 'analyze_clusters', args: { element: 'Ta', r_min: 2.65, r_max: 3.33 } }] },
            { text: 'In this one configuration, analyze_clusters reads the Ta4 clusters, not the TaSe6 octahedra.\nOutcome: partial — no distortion index.' }
        ], { run: demoRun, offered: shellTools.filter((tool) => tool.name !== 'polyhedral_distortion') });
        expect(done.gates).toEqual([{ round: 0, id: 'polyhedra', missing: ['analyze_clusters'] }]);
        expect(requests[1].at(-1).content).toMatch(/concludes a distortion of the coordination polyhedra, but this question has not run analyze_clusters\./);
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

// The user's interview answers of 2026-10-24 (local-symmetry), each where a
// check in code catches the wrong answer.
describe('local-symmetry: the user\'s interview answers (2026-10-24)', () => {
    const ladder = { calls: [{ name: 'detect_symmetry', args: { tolerance_A: 0.2 } }, { name: 'detect_symmetry', args: { tolerance_A: 0.02 } }] };
    const rulesOf = (done) => done.audit.map((finding) => finding.rule);

    it('answer 4a: a group–subgroup relation is read against the run\'s groups — the Demo\'s P-42₁m in F-43m: index 12, neither t nor k', async () => {
        const right = await replay([ladder, {
            text: 'From F-43m at 0.2 Å to P-42₁m at 0.02 Å the symmetry lowers: P-42₁m is a subgroup of F-43m of index 12, neither t nor k '
                + '(the F centring and the three-fold axes are lost).\nOutcome: achieved — both tolerances run.'
        }]);
        expect(right.toolMessages.map((message) => [message.data.space_group, message.data.n_ops])).toEqual([['F-43m', 96], ['P-42_1m', 8]]);
        expect(rulesOf(right.done)).toEqual([]);
        const wrong = await replay([ladder, { text: 'P-42₁m is a maximal t-subgroup of F-43m, of index 3.\nOutcome: achieved — x.' }]);
        expect(rulesOf(wrong.done)).toEqual(['subgroup']);
        expect(wrong.done.audit[0].text).toMatch(/^The answer has P-42_1m as a subgroup of F-43m, but its index is 96\/8 = 12 operations of the cell, not 3: /);
        expect(wrong.done.outcome.status).toBe('not verified');
    });

    it('answer 4b: a cell doubling the box cannot hold — the Demo\'s 10 × 10 × 10 box holds one, so only the superstructure test speaks there', async () => {
        const demoRun = { rValueFile, plotFiles, structure };
        expect(buildRunContext(demoRun).structure.supercell).toEqual([10, 10, 10]);
        const { done } = await replay([
            { calls: [{ name: 'analyze_modulations', args: {} }] },
            { text: 'In this configuration the Ta order doubles the cell along c.\nOutcome: achieved — x.' }
        ], { run: demoRun });
        expect(rulesOf(done)).toEqual(['superstructure']);
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
        expect(done.audit[0].text).toMatch(/Ga: peak -0\.887σ at [\d.]+× isotropic, map 0\.561σ — nothing to report/);
    });
});

// The user's interview answers of 2026-10-24 (displacement-anisotropy): the
// tiers on the Demo, and the disorder rule on a dev case (the Demo has no split
// or flat-topped site), each through the turn as the app runs it.
describe('displacement-anisotropy: the user\'s interview answers (2026-10-24)', () => {
    const rulesOf = (done) => done.audit.map((finding) => finding.rule);
    // One question on a dev case's run, as the harness opens it (never a held-out case).
    const devTurn = async (id, replies) => {
        const copilot = await loadRun(buildCase(DEV_CASES.find((item) => item.id === id)));
        const { callModel, requests } = scriptedModel([...replies]);
        const events = [];
        for await (const event of runTurn({ question: 'Q', run: copilot.props, model: callModel, tools: [...buildRunTools(copilot.props), ...copilot.tools] })) {
            events.push(event);
        }
        const toolMessages = requests.at(-1).filter((message) => message.role === 'tool').map((message) => JSON.parse(message.content));
        return { done: events.at(-1), toolMessages };
    };

    it('answers 3 and 4: Ta site 17 is a finding at two resolutions, Ta site 25 only a hint, and a preferred direction needs a finding', async () => {
        const both = { calls: [17, 25].map((ref) => ({ name: 'analyze_displacement_directions', args: { site_ref: ref } })) };
        const right = await replay([both, {
            text: 'Site 17 (Ta) has a preferred displacement direction: 3.15σ at 1.95× isotropic, toward {{r1.peak.direction_uvw}} ± 20.9° (cell radius), '
                + 'and it holds at ν = 1 (3.83σ). Site 25 (Ta) is only a hint: 4.47σ at ν = 2, but 2.33σ at ν = 1.\nOutcome: achieved — tiers.'
        }]);
        expect(right.toolMessages.map((message) => [message.data.site, message.data.peak.tier, message.data.neighbour?.holds])).toEqual([[17, 'finding', true], [25, 'hint', false]]);
        expect(rulesOf(right.done)).toEqual([]);
        // The hint alone called preferred: flagged; called a hint: not.
        const one = { calls: [{ name: 'analyze_displacement_directions', args: { site_ref: 25 } }] };
        const wrong = await replay([one, { text: 'Site 25 (Ta) has a preferred displacement direction at 4.47σ.\nOutcome: achieved — x.' }]);
        expect(rulesOf(wrong.done)).toEqual(['direction']);
        expect(wrong.done.audit[0].text).toMatch(/site 25: peak 4\.47σ at 2\.23× isotropic, map [\d.]+σ — a hint, not held at the neighbouring resolution/);
        const hint = await replay([one, { text: 'Site 25 (Ta) shows a hint of a preferred direction, not a finding: it does not hold at ν = 1.\nOutcome: partial — a hint.' }]);
        expect(rulesOf(hint.done)).toEqual([]);
    });

    it('answer 8 (firm rule): a split site called disorder before analyze_modulations ran for its element is flagged; after it, the order is read', async () => {
        // dev-zone-boundary's Ti: split (list_sites), and ordered at (0, 0, ½).
        const wrong = await devTurn('dev-zone-boundary', [
            { calls: [{ name: 'list_sites', args: { element: 'Ti' } }] },
            { text: 'Site 2 (Ti) shows positional disorder along c.\nOutcome: achieved — split.' }
        ]);
        expect(wrong.toolMessages[0].data.sites[0]).toMatchObject({ ref: 2 });
        expect(rulesOf(wrong.done)).toEqual(['disorder']);
        expect(wrong.done.outcome.status).toBe('not verified');
        const right = await devTurn('dev-zone-boundary', [
            { calls: [{ name: 'list_sites', args: { element: 'Ti' } }, { name: 'analyze_modulations', args: { element: 'Ti' } }] },
            { text: 'Site 2 (Ti) is split, but analyze_modulations finds Ti ordered at q = (0, 0, ½): ordered antiphase displacements doubling the cell '
                + 'along c, not disorder.\nOutcome: achieved — order.' }
        ]);
        expect(rulesOf(right.done)).toEqual([]);
    });
});

describe('local-order: the user\'s interview answers (2026-10-24)', () => {
    const demoRun = { rValueFile, plotFiles, structure };
    const orderTools = [...buildRunTools(demoRun), ...hostTools];
    const rulesOf = (done) => done.audit.map((finding) => finding.rule);
    const modulations = { calls: [{ name: 'analyze_modulations', args: {} }] };

    it('answer 6: the Demo\'s in-phase correlation is neighbouring cells moving together (correlated motion), not order — and never "no correlation"', async () => {
        const right = await replay([modulations, {
            text: 'Ga, Se and Ta displacements are correlated in phase between neighbouring cells: neighbouring cells moving together (correlated motion), '
                + 'not order, with no ordered zone-boundary point.\nOutcome: achieved — correlated motion.'
        }], { run: demoRun, offered: orderTools });
        const inPhase = right.toolMessages[0].checks.find((check) => /correlated in phase/.test(check.text));
        expect(inPhase.text).toMatch(/— neighbouring cells moving together \(correlated motion; as RMC's long waves and acoustic-like thermal motion give\), not order\.$/);
        expect(rulesOf(right.done)).toEqual([]);
        for (const [text, rule] of [
            ['Ga, Se and Ta show short-range order between cells.', 'correlated-motion'],
            ['The Ga, Se and Ta displacements are uncorrelated between cells.', 'no-correlation']
        ]) {
            const wrong = await replay([modulations, { text: `${text}\nOutcome: achieved — x.` }], { run: demoRun, offered: orderTools });
            expect({ text, rules: rulesOf(wrong.done) }).toEqual({ text, rules: [rule] });
            expect(wrong.done.outcome.status).toBe('not verified');
        }
    });

    // Runs the real modulation, direction and cluster engines on the 52,000-atom Demo three times, and the polyhedra engine once (over 5 s on a loaded CI runner).
    it('answer 5: "no order" needs polyhedral_distortion on the Demo\'s GaSe₄ only where the run offers it; its diffuse verdict is analyze_modulations\' own', async () => {
        const list = { calls: [
            { name: 'analyze_modulations', args: {} },
            { name: 'analyze_displacement_directions', args: { site_ref: 1 } },
            { name: 'analyze_clusters', args: { element: 'Ta', r_min: 2.65, r_max: 3.33 } }
        ] };
        const answer = { text: 'There is no local order beyond thermal motion in the tests that ran: in phase only, the Ga, Se and Ta correlation is '
            + 'correlated motion.\nOutcome: partial — the tests that ran.' };
        // Without polyhedral_distortion: the list as before, met.
        const withoutTools = orderTools.filter((tool) => tool.name !== 'polyhedral_distortion');
        const without = await replay([list, answer], { run: demoRun, offered: withoutTools });
        expect(without.done.gates).toEqual([]);
        expect(rulesOf(without.done)).toEqual([]);
        // With it offered: Ga's tetrahedra (no finite cluster) are asked for, and a call centred on Ga meets the list.
        const polyhedra = {
            name: 'polyhedral_distortion',
            description: 'Distortion of the coordination polyhedra of a centre element (a stand-in for the tool).',
            parameters: { type: 'object', properties: { centre: { type: 'string' }, ligand: { type: 'string' } }, required: ['centre', 'ligand'] },
            run: async () => ({ data: { measured: true }, checks: [] })
        };
        const asked = await replay([list, answer, { calls: [{ name: 'polyhedral_distortion', args: { centre: 'Ga', ligand: 'Se' } }] }, answer],
            { run: demoRun, offered: [...withoutTools, polyhedra] });
        expect(asked.done.gates).toEqual([{ round: 1, id: 'no-order', missing: ['polyhedral_distortion for Ga (centre Ga: its polyhedra form no finite cluster)'] }]);
        expect(rulesOf(asked.done)).toEqual([]);
        // The real tool, as the Demo offers it: the same gate, met by a call centred on Ga in its whole shell.
        const real = await replay([list, answer, { calls: [{ name: 'polyhedral_distortion', args: { centre: 'Ga', ligand: 'Se', r_min: 2.18, r_max: 2.79 } }] }, answer],
            { run: demoRun, offered: orderTools });
        expect(real.done.gates).toEqual([{ round: 1, id: 'no-order', missing: ['polyhedral_distortion for Ga (centre Ga: its polyhedra form no finite cluster)'] }]);
        expect(rulesOf(real.done)).toEqual([]);
    }, 30_000);
});

describe('the real model\'s correct answers (gemma4:26b, 2026-10-11)', () => {
    // The checks on an answer through the verification's interface: the
    // tools' facts recorded, the run context given, no claim offered.
    const auditAnswer = (text, { facts = [], context = null } = {}) => {
        const verification = createVerification({ context, offered: [] });
        facts.forEach((item) => verification.record({ name: 'analysis', ...toolResult({ facts: item }) }, () => ({ data: null })));
        return verification.judge(text).done.audit;
    };

    it('raise no check', () => {
        // The run context the answers were given on: still improving; F-43m over P-42₁m.
        const context = {
            convergence: { watchdog_status: 'improving' },
            symmetry: { space_group: 'F-43m', point_group: '-43m', ladder: [{ sg: 'P-42_1m', pg: '-42m', n_ops: 8 }, { sg: 'F-43m', pg: '-43m', n_ops: 96 }] }
        };
        expect(auditAnswer('The run has not converged; the dashboard\'s watchdog status is improving. The ln(χ²) of the last fit '
            + 'term (X_ray_(R)1) has decreased from 0.0953 to −8.87 over the course of the run, and the change over the last 10% '
            + 'of the run is −0.11, indicating that the fit is still actively descending.', { context })).toEqual([]);
        expect(auditAnswer('Based on the displacement analysis, the Ga atoms do not exhibit a significant preferred displacement '
            + 'direction. The pooled analysis for all Ga atoms shows no significant peak direction (significance of −0.887σ).', {
            facts: [{ direction: { subject: 'Ga', significant: false, peak_sigma: -0.887, map_sigma: 0.561 } }]
        })).toEqual([]);
        expect(auditAnswer('At a 0.2 Å tolerance, the space group is F-43m (No. 216). At a 0.02 Å tolerance, the space group is '
            + 'P-42_1m. The symmetry breaking between these two tolerances is the loss of the F-43m symmetry. Note that the '
            + 'transition from F-43m to P-42_1m occurs at 0.031 Å, which is the noise_floor_A for this run, meaning this change '
            + 'is within the sampling noise.', { context })).toEqual([]);
    });
});
