// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// What the model sees, 1 of 4: every system prompt the Copilot can build
// (copilotSystemPrompt in src/llm/prompts/copilot.js), one file per variant
// in __snapshots__/system-prompt/. The skills' index is the bundled
// .claude/skills (copilotSkills.js), which the validation harness reads from
// disk (validation/run.js loadSkills) into the same text. A change to these
// files is a change to what the model is told: call it out in the PR.

import { describe, expect, it } from 'vitest';
import { buildCopilotMessages, COPILOT_PROMPT, COPILOT_PROMPT_NO_TOOLS, COPILOT_PROMPT_VIEW } from '../../llm/prompts/copilot';
import { parseSkillFiles, skillIndex } from '../../llm/copilot/skills';
import { loadSkills } from '../../../validation/run';
import { bundledSkills, PROMPT_VARIANTS, promptVariants, renderPrompt, variantFlags } from './helpers';

const index = skillIndex(bundledSkills());
const variants = promptVariants(index);

describe('the system prompt the model receives', () => {
    it('reads the same skills\' index in the app and in the validation harness', () => {
        const disk = loadSkills();
        expect(skillIndex(parseSkillFiles(disk.files, { available: disk.references }))).toBe(index);
    });

    it.each(PROMPT_VARIANTS)('variant %s', async (name) => {
        await expect(renderPrompt(name, variants.get(name))).toMatchFileSnapshot(`__snapshots__/system-prompt/${name}.txt`);
    });

    it('is what buildCopilotMessages puts first, and the named prompts are variants', () => {
        for (const name of PROMPT_VARIANTS) {
            const flags = variantFlags(name);
            const [system] = buildCopilotMessages({ context: {}, userText: 'q', ...flags, skills: flags.skills ? index : '' });
            expect(system).toEqual({ role: 'system', content: variants.get(name) });
        }
        expect(COPILOT_PROMPT_NO_TOOLS).toBe(variants.get('no-tools'));
        expect(COPILOT_PROMPT).toBe(variants.get('tools'));
        expect(COPILOT_PROMPT_VIEW).toBe(variants.get('tools-view'));
        // With tools off nothing else changes the prompt.
        expect(buildCopilotMessages({ context: {}, userText: 'q', toolsEnabled: false, viewTools: true, askUser: true, skills: index })[0].content)
            .toBe(variants.get('no-tools'));
    });
});
