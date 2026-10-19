// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// What the model sees, at the turn module's interface (src/llm/copilot/turn.js,
// runTurn): the first request of a dock-shaped spec and of a harness-shaped
// spec is the request in the dock's and the harness's snapshot files. The turn
// is the one assembly behind both adapters, so each adapter adds only what is
// its own: the dock its run props, its tools (the run tools over its cards,
// highlight_cards, the host's engine and view tools), its question card and
// the user's settings through providerModel; the harness a case opened
// through the run module, its compute-only tools, the simulated user and the
// skills (arm B). The files are read, never written here; each file's two
// header lines name its own capture point and are not compared.

import { readFileSync } from 'node:fs';
import { dirname, join } from 'node:path';
import { fileURLToPath } from 'node:url';
import { beforeAll, describe, expect, it, vi } from 'vitest';
import { buildRunTools } from '../../llm/copilot/runTools';
import { skillIndex } from '../../llm/copilot/skills';
import { runTurn } from '../../llm/copilot/turn';
import { buildFocusTools } from '../../llm/copilot/viewTools';
import { providerModel } from '../../llm/provider/model';
import { DEFAULT_SETTINGS } from '../../llm/settings';
import { copilotSkills } from '../../copilotSkills';
import { buildCopilotCards, COPILOT_PAGES } from '../../copilotCards';
import { buildViewTools } from '../../copilotViewTools';
import { buildCase, DEV_CASES } from '../../../validation/cases';
import { liveQuestion, skippingUser } from '../../../validation/harness';
import { loadSkills } from '../../../validation/run';
import { loadRun } from '../../../validation/runFolder';
import { bundledSkills, jsonFile, loadDemoAssistantRun, promptVariants, recordingModel, renderRequest } from './helpers';

const provider = vi.hoisted(() => ({ current: null }));
vi.mock('../../llm/provider/contextGuard', async (importOriginal) => ({
    ...(await importOriginal()),
    streamCopilotChat: (request) => provider.current.callModel(request)
}));

const HERE = dirname(fileURLToPath(import.meta.url));
const snapshot = (name) => readFileSync(join(HERE, '__snapshots__', name), 'utf8');
// A rendered request without its two header lines (the title and the capture point).
const body = (text) => text.split('\n').slice(2).join('\n');
const variants = promptVariants(skillIndex(bundledSkills()));
const render = (request, options) => body(renderRequest(request, { title: '', port: '', variants, ...options }));

const run = async (spec) => {
    for await (const event of runTurn(spec)) void event;
};

describe('the turn\'s first request is the dock\'s and the harness\'s', () => {
    let demo;
    let built;
    let caseRun;

    beforeAll(async () => {
        demo = await loadDemoAssistantRun();
        built = buildCase(DEV_CASES.find((item) => item.id === 'dev-offcentre'));
        caseRun = await loadRun(built);
    });

    it('a dock-shaped spec (the Demo run): dock/first-request.txt and dock/tools.json', async () => {
        const { run: props, copilot } = demo;
        const cards = buildCopilotCards({ plotFiles: props.plotFiles, rValueFile: props.rValueFile, symmetry: props.symmetry });
        const hostTools = [...copilot.tools, ...buildViewTools({ cards, pages: COPILOT_PAGES, activePage: () => 'dashboard' })];
        const recorder = recordingModel();
        provider.current = recorder;
        await run({
            question: 'Has this run converged?',
            run: props,
            model: providerModel({ ...DEFAULT_SETTINGS, model: 'snapshot-model' }, { session: {} }),
            tools: [...buildRunTools(props, { cards }), ...buildFocusTools({ cards, pages: COPILOT_PAGES }), ...hostTools],
            askUser: async () => ({ source: 'skipped' }),
            skills: copilotSkills
        });
        const [first] = recorder.requests;
        expect(jsonFile(first.tools)).toBe(snapshot('dock/tools.json'));
        expect(render(first, { toolsFile: 'dock/tools.json' })).toBe(body(snapshot('dock/first-request.txt')));
    });

    it('a harness-shaped spec (dev-offcentre): harness/first-request.arm-B.txt, arm-C.txt and tools.arm-B.json', async () => {
        const firstRequest = async (skills) => {
            const recorder = recordingModel();
            await run({
                question: liveQuestion(built),
                run: caseRun.props,
                model: recorder.port,
                tools: [...buildRunTools(caseRun.props, { cards: [] }), ...caseRun.tools],
                askUser: skippingUser([]),
                skills
            });
            return recorder.requests[0];
        };
        const b = await firstRequest(loadSkills());
        const c = await firstRequest(null);
        expect(jsonFile(b.tools)).toBe(snapshot('harness/tools.arm-B.json'));
        expect(render(b, { toolsFile: 'harness/tools.arm-B.json' })).toBe(body(snapshot('harness/first-request.arm-B.txt')));
        expect(render(c, { toolsFile: 'harness/tools.arm-B.json, without read_skill', base: b, baseName: 'first-request.arm-B.txt' }))
            .toBe(body(snapshot('harness/first-request.arm-C.txt')));
    });
});
