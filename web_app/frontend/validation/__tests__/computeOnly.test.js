// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// The live arms offer compute-only tools: they read their tools through the
// one list of tool effects (src/llm/copilot/effects.js, liveArmTools), so a
// view or change tool never reaches them, and the turn they run gets no
// approve port, without which it offers no change tool. The protocol is
// unchanged: the arms' tools are the ones they had before the list (the
// model-visible snapshots, harness/tools.arm-*.json, pin the bytes). A compute
// tool added after the 0.8.0 live gate (protocol 2026-11-01.1) is left out by
// name (LIVE_ARM_EXCLUDED: check_ripples) until a protocol bump gives it to
// the arms. Dev cases only; the held-out cases are never built here.

import { describe, expect, it, vi } from 'vitest';
import { buildRunTools } from '../../src/llm/copilot/runTools';
import { COPILOT_TOOL_EFFECTS, effectOf } from '../../src/llm/copilot/effects';
import { toRequestTools } from '../../src/llm/copilot/tools';
import { buildCase, DEV_CASES } from '../cases';
import { LIVE_ARM_EFFECTS, LIVE_ARM_EXCLUDED, liveArmTools, PROTOCOL_VERSION, runLiveCase } from '../harness';
import { loadRun } from '../runFolder';
import { openRunFolder } from '../../src/openRunFolder';
import { demoEntries } from '../../src/__tests__/modelVisible/helpers';

const MODEL = { name: 'scripted', label: 'scripted test model', model: 'scripted', baseUrl: 'none', temperature: 0 };
const answer = 'No order.\n\nOutcome: achieved — read it.\n\nVERDICT {"features": [], "space_group": "P1", "chiral": true, "centrosymmetric": false, "polar": true}';
const scripted = () => async function* call() {
    yield { content: answer };
    yield { finishReason: 'stop' };
};

// A change tool as the dock offers one, and a view tool: neither may reach an arm.
const changeTool = () => ({
    name: 'set_bond_geometry', description: 'Set the Bond Geometry inputs.', parameters: { type: 'object', properties: {} },
    run: vi.fn(async () => ({ data: { applied: true }, checks: [] }))
});
const viewTool = () => ({ name: 'current_view', view: true, description: 'What is on screen.', parameters: { type: 'object', properties: {} }, run: vi.fn() });

describe('the live arms offer compute-only tools (the effects list)', () => {
    it('LIVE_ARM_EFFECTS is compute only', () => {
        expect(LIVE_ARM_EFFECTS).toEqual(['compute']);
    });

    it('on every dev case: the tools they had before the list, every one listed as compute', async () => {
        for (const item of DEV_CASES) {
            const run = await loadRun(buildCase(item));
            const tools = liveArmTools(run);
            // Each tool as the model reads it: its name, description and parameters.
            const before = [...buildRunTools(run.props, { cards: [] }), ...run.tools]
                .filter((tool) => !LIVE_ARM_EXCLUDED.includes(tool.name));
            expect(toRequestTools(tools), item.id).toEqual(toRequestTools(before));
            for (const tool of tools) {
                expect(COPILOT_TOOL_EFFECTS[tool.name], `${item.id}: ${tool.name}`).toBe('compute');
            }
        }
    // Every dev case is built and opened: well within this on a loaded machine.
    }, 60000);

    it('leave check_ripples out by name: the dock has it, the arms of protocol 2026-11-01.1 do not', async () => {
        expect(LIVE_ARM_EXCLUDED).toEqual(['check_ripples']);
        expect(PROTOCOL_VERSION).toBe('2026-11-01.1');
        expect(COPILOT_TOOL_EFFECTS.check_ripples).toBe('compute');
        // The Demo has fit files, so its run carries the tool; the arms drop it.
        const demo = await openRunFolder(demoEntries());
        expect(demo.tools.map((tool) => tool.name)).toContain('check_ripples');
        expect(liveArmTools(demo).map((tool) => tool.name)).not.toContain('check_ripples');
        for (const item of DEV_CASES) {
            const run = await loadRun(buildCase(item));
            expect(liveArmTools(run).map((tool) => tool.name), item.id).not.toContain('check_ripples');
        }
    }, 60000);

    it('drop a view or change tool that a run would carry', async () => {
        const run = await loadRun(buildCase(DEV_CASES[0]));
        const tools = liveArmTools({ ...run, tools: [...run.tools, changeTool(), viewTool()] });
        expect(tools.map((tool) => tool.name)).toEqual(liveArmTools(run).map((tool) => tool.name));
        expect(tools.some((tool) => effectOf(tool) !== 'compute')).toBe(false);
    });

    it('a run offers no view or change tool: the record\'s tools are compute only, and no change ran', async () => {
        const record = await runLiveCase({ built: buildCase(DEV_CASES[0]), arm: 'B', model: MODEL, port: scripted() });
        expect(record.error).toBeNull();
        expect(record.tools_offered.length).toBeGreaterThan(5);
        expect(record.tools_offered.map((name) => effectOf({ name }))).toEqual(record.tools_offered.map(() => 'compute'));
        expect(record.tools_offered).not.toContain('set_bond_geometry');
    });
});
