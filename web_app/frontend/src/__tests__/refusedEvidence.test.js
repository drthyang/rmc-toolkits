// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// A refused analysis backs no claim (the user's decision A5, 2026-10-24),
// against the real tools on the Demo run: the evidence gate reads a window cut
// from the refusal's data (`refused`, the field compute_bond_angles and
// polyhedral_distortion declare in their dataShape) and offers both fixes —
// the whole shell, or the user's own window confirmed — and says any other
// refusal only as a refusal. The engines run in-thread (pcaKdeWorker's
// inProcessEngine, as src/openRunFolder.js runs them).

import { beforeAll, describe, expect, it } from 'vitest';
import { declaredKeys } from '../llm/copilot/dataShape';
import { buildRunTools } from '../llm/copilot/runTools';
import { executeToolCall } from '../llm/copilot/tools';
import { createVerification } from '../llm/copilot/verification';
import { buildCopilotTools } from '../copilotTools';
import { inProcessEngine } from '../workers/pcaKdeWorker';
import { loadDemoAssistantRun } from './modelVisible/helpers';

const CUT = { end1: 'Se', apex: 'Ta', end2: 'Se', r12_min: 2.3, r12_max: 2.7 };
const WHOLE = { end1: 'Se', apex: 'Ta', end2: 'Se', r12_min: 2.21, r12_max: 3.23 };
const EMPTY = { end1: 'Se', apex: 'Ta', end2: 'Se', r12_min: 0.5, r12_max: 1 };
const POLYHEDRA_CUT = { centre: 'Ta', ligand: 'Se', r_min: 2.3, r_max: 2.7 };
const WINDOW_FIXES = 'Run it again on the whole shell, or, if the user wants that window, ask them (ask_user) and call again with window_confirmed: true;';

describe('a refused analysis backs no claim, on the real tools (the Demo run)', () => {
    let tools;
    const results = {};
    const call = async (name, args) => {
        const tool = tools.find((candidate) => candidate.name === name);
        return { name, interactive: false, reading: false, view: false, ...(await executeToolCall(tool, JSON.stringify(args), { timeoutMs: 120000 })) };
    };
    // The send-back an answer gets after `steps`, while tools can still run.
    const sendBack = (text, steps) => {
        const verification = createVerification();
        steps.forEach((step) => verification.record(step, () => ({ data: step.data })));
        return verification.judge(text, { canSendBack: true }).sendBack ?? null;
    };

    beforeAll(async () => {
        const { run, rmc6fText } = await loadDemoAssistantRun();
        tools = [
            ...buildRunTools(run, { cards: [] }),
            ...buildCopilotTools({
                request: inProcessEngine(rmc6fText), structure: run.structure, elements: Object.keys(run.structure.elementCounts || {}).sort(),
                plotFiles: run.plotFiles, runSettings: run.runSettings
            })
        ];
        results.cut = await call('compute_bond_angles', CUT);
        results.whole = await call('compute_bond_angles', WHOLE);
        results.empty = await call('compute_bond_angles', EMPTY);
        results.polyhedraCut = await call('polyhedral_distortion', POLYHEDRA_CUT);
    }, 120000);

    it('the window-cut refusals carry `refused` in their data; a window with no angles refuses without it', () => {
        expect(results.cut).toMatchObject({ resultStatus: 'refused', error: null });
        expect(typeof results.cut.data.refused).toBe('string');
        expect(results.polyhedraCut).toMatchObject({ resultStatus: 'refused', error: null });
        expect(typeof results.polyhedraCut.data.refused).toBe('string');
        expect(results.empty).toMatchObject({ resultStatus: 'refused', error: null });
        expect(results.empty.data.refused).toBeUndefined();
        expect(results.whole).toMatchObject({ resultStatus: 'ok' });
    });

    it('the tools whose data can say `refused` are exactly the ones with window_confirmed: the window fixes are offered only where they exist', () => {
        const declaring = tools.filter((tool) => tool.dataShape && declaredKeys(tool.dataShape).includes('refused')).map((tool) => tool.name);
        const confirmable = tools.filter((tool) => tool.parameters?.properties?.window_confirmed).map((tool) => tool.name);
        expect(declaring).toEqual(['compute_bond_angles', 'polyhedral_distortion']);
        expect(confirmable).toEqual(declaring);
    });

    it('a coordination number on the cut window goes back with both fixes; on the whole shell it stands', () => {
        const sent = sendBack('Each Ta has CN 6.', [results.cut]);
        expect(sent).toMatchObject({ id: 'coordination', missing: ['find_gr_peaks'] });
        expect(sent.nudge).toContain('but compute_bond_angles refused its window (it cuts the coordination shell), and a refused call is no evidence.');
        expect(sent.nudge).toContain(WINDOW_FIXES);
        expect(sendBack('Each Ta has CN 6.', [results.cut, results.whole])).toBeNull();
    });

    it('distorted octahedra on polyhedral_distortion\'s cut window go back with both fixes', () => {
        const sent = sendBack('The TaSe6 octahedra are distorted.', [results.polyhedraCut]);
        expect(sent).toMatchObject({ id: 'polyhedra', missing: ['polyhedral_distortion'] });
        expect(sent.nudge).toContain('but polyhedral_distortion refused its window (it cuts the coordination shell)');
        expect(sent.nudge).toContain(WINDOW_FIXES);
    });

    it('a window with no angles is said as a refusal only: no window to confirm', () => {
        const sent = sendBack('Each Ta has CN 6.', [results.empty]);
        expect(sent.nudge).toContain('but compute_bond_angles refused, and a refused call is no evidence. Run it so that it answers,');
        expect(sent.nudge).not.toContain('window_confirmed');
    });
});
