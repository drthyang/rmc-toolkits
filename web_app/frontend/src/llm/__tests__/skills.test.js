// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

import { describe, expect, it, vi } from 'vitest';
import { runCopilotTurn } from '../copilot/agentLoop';
import {
    DRAFT_MARK, READ_SKILL, SKILL_RESULT_BUDGET, buildReadSkillTool, frontmatter, joinPath, parseSkill, parseSkillFiles, skillIndex
} from '../copilot/skills';
import { copilotSystemPrompt } from '../prompts/copilot';

const skillFile = (name, { draft = true, body = '' } = {}) => [
    '---',
    `name: ${name}`,
    'description: >-',
    `  The ${name} method, in one`,
    '  folded line.',
    '---',
    '',
    `# ${name}`,
    '',
    draft ? `**Status: ${DRAFT_MARK}.**` : '**Status: the user\'s method.**',
    '',
    body || '## Steps\n\n1. Do it.',
    '',
    '## References',
    '',
    `- [notes](notes.md) — what the ${name} notes cover.`
].join('\n');

const files = {
    '.claude/skills/shells/SKILL.md': skillFile('shells', { body: '## Steps\n\n1. Find the shell (to confirm with the user).' }),
    '.claude/skills/fit/SKILL.md': skillFile('fit', { draft: false })
};
const available = ['.claude/skills/shells/notes.md', '.claude/skills/fit/notes.md'];

describe('skill files', () => {
    it('read the frontmatter, folded blocks included', () => {
        expect(frontmatter('name: a\ndescription: >-\n  one\n  two\nother: "q"')).toEqual({ name: 'a', description: 'one two', other: 'q' });
        expect(joinPath('.claude/skills/a', '../b/x.md')).toBe('.claude/skills/b/x.md');
    });

    it('parse a skill: name, description, draft status, references beside it', () => {
        const skill = parseSkill(files['.claude/skills/shells/SKILL.md'], '.claude/skills/shells');
        expect(skill).toMatchObject({ name: 'shells', description: 'The shells method, in one folded line.', draft: true });
        expect(skill.references).toEqual([{ name: 'notes', path: '.claude/skills/shells/notes.md', description: 'what the shells notes cover.' }]);
        expect(() => parseSkill('# no frontmatter', 'd')).toThrow('has no frontmatter');
        expect(() => parseSkill('---\nname: x\n---\n', 'd')).toThrow('needs a name and a description');
    });

    it('parse every SKILL.md, sorted, and refuse a misnamed one or an unbundled reference', () => {
        expect(parseSkillFiles(files, { available }).map((skill) => skill.name)).toEqual(['fit', 'shells']);
        expect(() => parseSkillFiles({ '.claude/skills/x/SKILL.md': skillFile('y') })).toThrow('the name must be its folder\'s');
        expect(() => parseSkillFiles(files, { available: [] })).toThrow(/cites \.claude\/skills\/(fit|shells)\/notes\.md, which is not bundled/);
    });

    it('list only names and descriptions in the system prompt', () => {
        const skills = parseSkillFiles(files, { available });
        const index = skillIndex(skills);
        expect(index).toBe('- fit: The fit method, in one folded line.\n- shells: The shells method, in one folded line.');
        const prompt = copilotSystemPrompt({ skills: index });
        expect(prompt).toContain('read it with read_skill');
        expect(prompt).toContain(index);
        expect(prompt).not.toContain('Find the shell');
        expect(copilotSystemPrompt({ toolsEnabled: false, skills: index })).not.toContain('read_skill');
    });
});

describe('read_skill', () => {
    const skills = parseSkillFiles(files, { available });
    const load = vi.fn(async (path) => `# Notes\n\nFrom ${path}.`);
    const tool = buildReadSkillTool(skills, { load });

    it('is a kept, long-budget tool over the skills\' names', () => {
        expect(tool).toMatchObject({ name: READ_SKILL, keep: true, resultBudget: SKILL_RESULT_BUDGET });
        expect(tool.parameters.properties.name.enum).toEqual(['fit', 'shells']);
        expect(buildReadSkillTool([])).toBeNull();
    });

    it('returns a skill\'s body with its references listed, and a draft note', async () => {
        const { data, checks } = await tool.run({ name: 'shells' });
        expect(data.skill).toBe('shells');
        expect(data.text).toContain('Find the shell (to confirm with the user).');
        expect(data.text).toMatch(/References \(read one with read_skill and `reference`\):\n- notes: what the shells notes cover\.$/);
        expect(checks[0].text).toMatch(/draft awaiting the user's interview/);
        expect((await tool.run({ name: 'fit' })).checks).toEqual([]);
    });

    it('loads a reference on demand, and names the choices on a wrong name', async () => {
        const { data } = await tool.run({ name: 'shells', reference: 'notes' });
        expect(load).toHaveBeenCalledWith('.claude/skills/shells/notes.md');
        expect(data).toEqual({ skill: 'shells', reference: 'notes', text: '# Notes\n\nFrom .claude/skills/shells/notes.md.' });
        await expect(tool.run({ name: 'shells', reference: 'other' })).rejects.toThrow('shells has no reference "other"; it lists notes');
        await expect(tool.run({ name: 'nope' })).rejects.toThrow('no skill named "nope": the skills are fit, shells');
    });

    it('keeps a skill whole for the question: its own budget, never compacted', async () => {
        const long = { '.claude/skills/long/SKILL.md': skillFile('long', { body: `## Steps\n\n${'step. '.repeat(1000)}` }) };
        const reader = buildReadSkillTool(parseSkillFiles(long), {});
        const big = {
            name: 'big',
            description: 'Big.',
            parameters: { type: 'object', properties: { n: { type: 'number' } } },
            run: async ({ n }) => ({ data: { n, rows: Array.from({ length: 400 }, (_, i) => i) }, checks: [{ level: 'ok', text: `big ${n}` }] })
        };
        const replies = [
            [{ name: READ_SKILL, args: { name: 'long' } }],
            [{ name: 'big', args: { n: 1 } }],
            [{ name: 'big', args: { n: 2 } }],
            null
        ];
        const requests = [];
        const callModel = (request) => {
            requests.push(request.messages.map((message) => ({ ...message })));
            const calls = replies.shift();
            return (async function* stream() {
                if (!calls) { yield { content: 'Done.' }; return; }
                for (const [index, call] of calls.entries()) yield { toolCall: { index, id: `c${requests.length}`, name: call.name, arguments: JSON.stringify(call.args) } };
            }());
        };
        for await (const event of runCopilotTurn({
            callModel, messages: [{ role: 'system', content: 'S' }], tools: [reader, big], limits: { resultCharBudget: 1200, turnCharBudget: 1500 }
        })) void event;
        const toolMessages = requests.at(-1).filter((message) => message.role === 'tool');
        const skill = JSON.parse(toolMessages[0].content);
        // 6,000 characters of steps arrive whole, under its own budget, and stay so.
        expect(skill.data.text.length).toBeGreaterThan(6000);
        expect(skill).not.toHaveProperty('compacted');
        // The analysis results are fitted to the turn's budget around it.
        expect(JSON.parse(toolMessages[1].content).compacted).toMatch(/shortened/);
    });
});
