// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// Every Copilot tool has its effect in one place: the list in
// llm/copilot/effects.js (compute / view / change). On the Demo run, every
// tool the builders make — the tools over the run props, highlight_cards,
// the host's engine, view and change tools as App passes them
// (copilotHostTools.js), ask_user and read_skill — is named in the list, and
// the list names no other tool; `view: true` is set exactly on the view tools
// (the loop's and the prompt's switch), no tool declares an effect of its
// own, and a change tool alone has an `undo`. The dock reads its tools
// through the list with every effect, and offers the change tool only with
// its approve port; the validation's arms read theirs compute-only
// (validation/__tests__/computeOnly.test.js), so no change tool reaches them.
// What each adapter sends is pinned in modelVisible/ (dock.test.jsx,
// harness.test.js).

import { beforeAll, describe, expect, it } from 'vitest';
import { buildAskUserTool } from '../llm/copilot/askUser';
import { COPILOT_TOOL_EFFECTS, DOCK_EFFECTS, effectOf, toolsWithEffects } from '../llm/copilot/effects';
import { buildRunTools } from '../llm/copilot/runTools';
import { buildReadSkillTool } from '../llm/copilot/skills';
import { describeTurn } from '../llm/copilot/turn';
import { buildFocusTools } from '../llm/copilot/viewTools';
import { createChangeHistory } from '../bondGeometryControl';
import { buildCopilotCards, COPILOT_PAGES } from '../copilotCards';
import { copilotHostTools } from '../copilotHostTools';
import { copilotSkills } from '../copilotSkills';
import { LIVE_ARM_EFFECTS } from '../../validation/harness';
import { bundledSkills, loadDemoAssistantRun } from './modelVisible/helpers';

describe('one list of tool effects (llm/copilot/effects.js)', () => {
    let demo;
    let cards;
    let hostTools;
    let all;

    beforeAll(async () => {
        demo = await loadDemoAssistantRun();
        const { run, copilot } = demo;
        cards = buildCopilotCards({ plotFiles: run.plotFiles, rValueFile: run.rValueFile, symmetry: run.symmetry });
        hostTools = copilotHostTools({
            engineTools: copilot.tools, elements: copilot.elements, cards, pages: COPILOT_PAGES,
            activePage: () => 'dashboard', history: createChangeHistory()
        });
        all = [
            ...buildRunTools(run, { cards }),
            ...buildFocusTools({ cards, pages: COPILOT_PAGES }),
            ...hostTools,
            buildAskUserTool(),
            buildReadSkillTool(bundledSkills(), { load: copilotSkills.load })
        ];
    });

    it('names every tool the builders make on the Demo, and no other', () => {
        const names = all.map((tool) => tool.name);
        expect(new Set(names).size).toBe(names.length);
        expect([...names].sort()).toEqual(Object.keys(COPILOT_TOOL_EFFECTS).sort());
    });

    it('view: true exactly on the view tools; no tool declares its own effect; only a change tool undoes', () => {
        for (const tool of all) {
            const effect = COPILOT_TOOL_EFFECTS[tool.name];
            expect(Boolean(tool.view), tool.name).toBe(effect === 'view');
            expect(tool.effect, tool.name).toBeUndefined();
            expect(typeof tool.undo === 'function', tool.name).toBe(effect === 'change');
            expect(effectOf(tool), tool.name).toBe(effect);
        }
    });

    it('the dock reads every effect: the change tool offered with its approve port, never without one', () => {
        const dockTools = toolsWithEffects([
            ...buildRunTools(demo.run, { cards }),
            ...buildFocusTools({ cards, pages: COPILOT_PAGES }),
            ...hostTools
        ], DOCK_EFFECTS);
        const offered = (approve) => describeTurn({
            run: demo.run, tools: dockTools, askUser: async () => ({ source: 'skipped' }), approve, skills: copilotSkills
        }).offered.tools.map((tool) => tool.name);
        const withPort = offered(async () => ({ approved: false, auto: false }));
        expect(withPort).toContain('set_bond_geometry');
        expect(new Set(withPort.map((name) => COPILOT_TOOL_EFFECTS[name]))).toEqual(new Set(['compute', 'view', 'change']));
        expect(offered(null)).toEqual(withPort.filter((name) => name !== 'set_bond_geometry'));
    });

    it('the arms read compute only: none of the host\'s view or change tools', () => {
        expect(LIVE_ARM_EFFECTS).toEqual(['compute']);
        const armTools = toolsWithEffects(hostTools, LIVE_ARM_EFFECTS).map((tool) => tool.name);
        expect(armTools).toEqual(demo.copilot.tools.map((tool) => tool.name));
        expect(armTools).not.toContain('set_bond_geometry');
        expect(armTools).not.toContain('orient_3d_view');
    });
});
