// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// What the model sees, 1 of 4: every system prompt the Copilot can build
// (copilotSystemPrompt in src/llm/prompts/copilot.js), one file per variant
// in __snapshots__/system-prompt/, those with method 'off' (the validation's
// arm C, llm copilot/method.js) named *.method-off. The skills' index is the bundled
// .claude/skills (copilotSkills.js), which the validation harness reads from
// disk (validation/run.js loadSkills) into the same text. A change to these
// files is a change to what the model is told: call it out in the PR.

import { describe, expect, it } from 'vitest';
import {
    buildCopilotMessages, COPILOT_PROMPT, COPILOT_PROMPT_NO_TOOLS, COPILOT_PROMPT_VIEW, EVIDENCE_RULES, METHOD_SENTENCES
} from '../../llm/prompts/copilot';
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
        // With tools off nothing else changes the prompt, the method included.
        expect(buildCopilotMessages({ context: {}, userText: 'q', toolsEnabled: false, viewTools: true, askUser: true, skills: index })[0].content)
            .toBe(variants.get('no-tools'));
        expect(buildCopilotMessages({ context: {}, userText: 'q', toolsEnabled: false, askUser: true, method: 'off' })[0].content)
            .toBe(variants.get('no-tools'));
    });

    it('method off (arm C): each variant is its full one without EVIDENCE_RULES and the method sentences, and never the skills\' index', () => {
        const off = PROMPT_VARIANTS.filter((name) => name.endsWith('.method-off'));
        expect(off).toEqual(['tools.method-off', 'tools-ask.method-off', 'tools-view.method-off', 'tools-view-ask.method-off']);
        for (const name of off) {
            const text = variants.get(name);
            expect(text).not.toContain(EVIDENCE_RULES);
            expect(text).not.toContain('Skills —');
            METHOD_SENTENCES.forEach((sentence) => expect(text).not.toContain(sentence));
            let full = variants.get(name.replace('.method-off', ''));
            [EVIDENCE_RULES, ...METHOD_SENTENCES].forEach((sentence) => { full = full.replace(` ${sentence}`, ''); });
            expect(text).toBe(full);
            // A skills' index passed with method off is not sent.
            const flags = variantFlags(name);
            expect(buildCopilotMessages({ context: {}, userText: 'q', ...flags, skills: index })[0].content).toBe(text);
        }
    });
});
