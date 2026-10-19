// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// What the model sees, 2-3 of 4, the validation harness's path: the first
// request runLiveCase (validation/harness.js) sends on a dev case, captured at
// its model port (`port`: a recording port here, which records the turn's
// intent — messages, tools, tool choice — with its signal; validation/run.js
// gives providerModel, which adds the model's base URL, id and temperature
// and sends through streamCopilotChat), with the simulated user that skips
// every question. No model runs.
//
//   harness/first-request.arm-B.txt   arm B (skills): fields, tool names, messages
//   harness/first-request.arm-C.txt   arm C (no skills), against arm B
//   harness/tools.arm-B.json          the `tools` field exactly as sent (arm B)
//
// A change to these files is a change to what the model sees: call it out in
// the PR. The files' header line still names the capture point by its name
// before the turn module (runLiveCase's callModel): PORT below is kept as it
// was so the files stay byte-identical.

import { describe, expect, it } from 'vitest';
import { buildCase, DEV_CASES } from '../../../validation/cases';
import { runLiveCase } from '../../../validation/harness';
import { loadSkills } from '../../../validation/run';
import { skillIndex } from '../../llm/copilot/skills';
import { bundledSkills, jsonFile, promptVariants, recordingModel, renderRequest } from './helpers';

const CASE = 'dev-offcentre';
const MODEL = { name: 'snapshot', label: 'recording fake model', model: 'snapshot', baseUrl: 'none', temperature: 0.2 };
const variants = promptVariants(skillIndex(bundledSkills()));
const PORT = 'runLiveCase\'s callModel(request) (validation/harness.js), which validation/run.js binds to '
    + 'streamCopilotChat({ baseUrl, model, temperature, apiKey, session, ...request })';

const firstRequest = async (arm) => {
    const model = recordingModel();
    const record = await runLiveCase({
        built: buildCase(DEV_CASES.find((item) => item.id === CASE)), arm, model: MODEL,
        port: model.port, skills: arm === 'B' ? loadSkills() : null
    });
    expect(record.error).toBeNull();
    return { request: model.requests[0], record };
};

describe(`the harness's first request (${CASE})`, () => {
    it('arm B: the request and its tools', async () => {
        const { request, record } = await firstRequest('B');
        expect(request.tools.map((tool) => tool.function.name)).toEqual(record.tools_offered);
        await expect(jsonFile(request.tools)).toMatchFileSnapshot('__snapshots__/harness/tools.arm-B.json');
        await expect(renderRequest(request, {
            title: `First request of the validation harness, arm B (skills), dev case ${CASE}`,
            port: PORT, variants, toolsFile: 'harness/tools.arm-B.json'
        })).toMatchFileSnapshot('__snapshots__/harness/first-request.arm-B.txt');
    });

    it('arm C: the same request without the skills', async () => {
        const [{ request: b }, { request: c }] = [await firstRequest('B'), await firstRequest('C')];
        // Arm C's tools are arm B's without read_skill, each definition unchanged.
        expect(c.tools).toEqual(b.tools.filter((tool) => tool.function.name !== 'read_skill'));
        await expect(renderRequest(c, {
            title: `First request of the validation harness, arm C (no skills), dev case ${CASE}`,
            port: PORT, variants, toolsFile: 'harness/tools.arm-B.json, without read_skill',
            base: b, baseName: 'first-request.arm-B.txt'
        })).toMatchFileSnapshot('__snapshots__/harness/first-request.arm-C.txt');
    });
});
