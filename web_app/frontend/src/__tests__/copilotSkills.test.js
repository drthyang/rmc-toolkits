// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// The repository's Copilot skills (.claude/skills, shared with Claude Code)
// as the app bundles them: each parses, is a draft awaiting the user's
// interview or an interviewed skill (its status line names the day;
// docs/SKILL_INTERVIEWS.md records the answers), has the sections a skill
// needs, and reaches the agent whole.

import { describe, expect, it } from 'vitest';
import { copilotSkills } from '../copilotSkills';
import { DRAFT_MARK, SKILL_RESULT_BUDGET, buildReadSkillTool, parseSkillFiles, skillIndex } from '../llm/copilot/skills';

const skills = parseSkillFiles(copilotSkills.files, { available: copilotSkills.references });
const reader = buildReadSkillTool(skills, { load: copilotSkills.load });
const INTERVIEWED = /\*\*Status: Interviewed \d{4}-\d{2}-\d{2}\.\*\*/;

describe('the Copilot\'s skills', () => {
    it('are the four drafted from the spot check and the review, and local-order from validation batch 1', () => {
        expect(skills.map((skill) => skill.name)).toEqual([
            'coordination-shells', 'displacement-anisotropy', 'local-order', 'local-symmetry', 'rmc-fit-assessment'
        ]);
    });

    it('are interviewed one by one (docs/SKILL_INTERVIEWS.md)', () => {
        expect(skills.filter((skill) => INTERVIEWED.test(skill.body)).map((skill) => skill.name)).toEqual(['coordination-shells', 'displacement-anisotropy', 'local-symmetry', 'rmc-fit-assessment']);
    });

    it.each(skills.map((skill) => [skill.name, skill]))('%s: a draft or an interviewed skill, with steps, a checklist and failure modes', (_, skill) => {
        const interviewed = INTERVIEWED.test(skill.body);
        expect(skill.draft).toBe(!interviewed);
        if (interviewed) {
            // The user's answers replaced every open choice.
            expect(skill.body).not.toMatch(/\(to confirm with the user\)/);
        } else {
            expect(skill.body).toContain(`**Status: ${DRAFT_MARK}.**`);
            expect(skill.body).toContain('## Judgment calls');
            expect(skill.body).toMatch(/\(to confirm with the user\)/);
        }
        for (const section of ['## Steps', '## Acceptance checklist', '## Known failure modes', '## Report', '## References']) {
            expect(skill.body).toContain(section);
        }
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
