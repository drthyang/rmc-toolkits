// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// The repository's Copilot skills (.claude/skills, shared with Claude Code)
// as the app bundles them: each parses, is a draft awaiting the user's
// interview, has the sections a skill needs, and reaches the agent whole.

import { describe, expect, it } from 'vitest';
import { copilotSkills } from '../copilotSkills';
import { DRAFT_MARK, SKILL_RESULT_BUDGET, buildReadSkillTool, parseSkillFiles, skillIndex } from '../llm/copilot/skills';

const skills = parseSkillFiles(copilotSkills.files, { available: copilotSkills.references });
const reader = buildReadSkillTool(skills, { load: copilotSkills.load });

describe('the Copilot\'s skills', () => {
    it('are the four drafted from the spot check and the review, and local-order from validation batch 1', () => {
        expect(skills.map((skill) => skill.name)).toEqual([
            'coordination-shells', 'displacement-anisotropy', 'local-order', 'local-symmetry', 'rmc-fit-assessment'
        ]);
    });

    it.each(skills.map((skill) => [skill.name, skill]))('%s: a draft with steps, judgment calls, a checklist and failure modes', (_, skill) => {
        expect(skill.draft).toBe(true);
        expect(skill.body).toContain(`**Status: ${DRAFT_MARK}.**`);
        for (const section of ['## Steps', '## Judgment calls', '## Acceptance checklist', '## Known failure modes', '## Report', '## References']) {
            expect(skill.body).toContain(section);
        }
        expect(skill.body).toMatch(/\(to confirm with the user\)/);
        expect(skill.body).toMatch(/- \[ \] /);
        expect(skill.body).toMatch(/\*Firm rule in\s+the app:\*/);
        // One line in the system prompt, not the text.
        expect(skill.description.length).toBeLessThan(600);
        expect(skill.references.length).toBeGreaterThan(0);
    });

    it.each(skills.map((skill) => [skill.name]))('%s: reaches the agent whole, with its references', async (name) => {
        const { data } = await reader.run({ name });
        expect(JSON.stringify({ checks: [], data }).length).toBeLessThan(SKILL_RESULT_BUDGET);
        const skill = skills.find((item) => item.name === name);
        for (const reference of skill.references) {
            const loaded = await reader.run({ name, reference: reference.name });
            expect(loaded.data.text.length).toBeGreaterThan(200);
            expect(JSON.stringify(loaded).length).toBeLessThan(SKILL_RESULT_BUDGET);
        }
    });

    it('cost the system prompt only their names and descriptions', () => {
        const index = skillIndex(skills);
        expect(index.split('\n')).toHaveLength(5);
        expect(index.length).toBeLessThan(2800);
    });
});
