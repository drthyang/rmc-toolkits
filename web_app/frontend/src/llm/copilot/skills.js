// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// The Copilot's skills: the methods it follows, written as Claude Code skills
// (`.claude/skills/<name>/SKILL.md` at the repository root, frontmatter `name`
// and `description`), so Claude Code and the in-app agent read the same text.
//
// Progressive disclosure: the system prompt lists each skill's name and
// description only (skillIndex); the agent reads a skill's body with the
// read_skill tool when a question calls for it, and a reference (listed in the
// skill's "## References" section, a file beside SKILL.md) the same way. The
// fixed prompt stays small, which a local model's context needs.
//
// The module never imports the files itself (its import boundary): the host
// bundles them and passes `{ files, load }` — the SKILL.md texts by repository
// path, and a loader for the reference files (copilotSkills.js).

export const READ_SKILL = 'read_skill';
// A skill is longer than an analysis result; it gets its own budget, and it is
// never compacted away later in the question (the agent follows it throughout).
export const SKILL_RESULT_BUDGET = 9000;
export const DRAFT_MARK = 'Draft: awaiting the user\'s interview';

/**
 * The frontmatter's top-level scalars: `key: value`, a quoted value, or a
 * folded / literal block (`>-`, `|`) of indented lines. Nested keys are
 * ignored; only the name and the description are needed.
 */
export const frontmatter = (text) => {
    const out = {};
    const lines = text.split(/\r?\n/);
    for (let i = 0; i < lines.length; i += 1) {
        const pair = /^([A-Za-z][\w-]*):\s*(.*)$/.exec(lines[i]);
        if (!pair) continue;
        const [, key, raw] = pair;
        if (/^[>|][-+]?$/.test(raw)) {
            const block = [];
            while (i + 1 < lines.length && (/^\s+\S/.test(lines[i + 1]) || lines[i + 1] === '')) {
                i += 1;
                block.push(lines[i].trim());
            }
            out[key] = raw.startsWith('>') ? block.join(' ').replace(/\s+/g, ' ').trim() : block.join('\n').trim();
        } else {
            out[key] = raw.replace(/^(["'])(.*)\1$/, '$2').trim();
        }
    }
    return out;
};

/** `dir` joined with a relative path, `..` resolved: a repository path. */
export const joinPath = (dir, relative) => {
    const out = [];
    for (const part of [...dir.split('/'), ...relative.split('/')]) {
        if (part === '' || part === '.') continue;
        if (part === '..') out.pop();
        else out.push(part);
    }
    return out.join('/');
};

// The "## References" section: `- [label](relative/path.md) — what it covers`.
const referencesOf = (body, dir) => {
    const section = /^## References\s*\n([\s\S]*?)(?=^## |(?![\s\S]))/m.exec(body);
    if (!section) return [];
    const out = [];
    for (const line of section[1].split('\n')) {
        const ref = /^\s*[-*]\s+\[([^\]]+)\]\(([^)]+)\)\s*(?:[—–-]+\s*)?(.*)$/.exec(line);
        if (!ref) continue;
        out.push({ name: ref[1].replace(/`/g, '').trim(), path: joinPath(dir, ref[2].trim()), description: ref[3].trim() });
    }
    return out;
};

/** One SKILL.md: `{ name, description, body, draft, references }`. Throws on a malformed file. */
export const parseSkill = (source, dir) => {
    const match = /^---\r?\n([\s\S]*?)\r?\n---\r?\n?/.exec(source);
    if (!match) throw new Error(`${dir}/SKILL.md has no frontmatter`);
    const meta = frontmatter(match[1]);
    if (!meta.name || !meta.description) throw new Error(`${dir}/SKILL.md needs a name and a description`);
    const body = source.slice(match[0].length).trim();
    return {
        name: meta.name,
        description: meta.description,
        body,
        draft: body.includes(DRAFT_MARK),
        references: referencesOf(body, dir)
    };
};

/**
 * Every skill in `files` (`{ '.claude/skills/<name>/SKILL.md': text }`),
 * sorted by name. A reference that is not among `available` paths (when
 * given) is an error: the in-app agent could not read it.
 */
export const parseSkillFiles = (files, { available = null } = {}) => Object.entries(files || {})
    .filter(([path]) => /^\.claude\/skills\/[^/]+\/SKILL\.md$/.test(path))
    .map(([path, source]) => {
        const dir = path.replace(/\/SKILL\.md$/, '');
        const skill = parseSkill(source, dir);
        if (`.claude/skills/${skill.name}` !== dir) throw new Error(`${dir}/SKILL.md is named "${skill.name}": the name must be its folder's`);
        if (available) {
            const missing = skill.references.find((reference) => !available.includes(reference.path));
            if (missing) throw new Error(`${dir}/SKILL.md cites ${missing.path}, which is not bundled`);
        }
        return skill;
    })
    .sort((a, b) => a.name.localeCompare(b.name));

/** The system prompt's list: one line per skill, its name and description. */
export const skillIndex = (skills) => (skills || []).map((skill) => `- ${skill.name}: ${skill.description}`).join('\n');

/**
 * read_skill(name, reference?): a skill's body (with its references listed),
 * or one reference's text. `load(path)` reads a bundled reference file.
 */
export const buildReadSkillTool = (skills, { load = null } = {}) => {
    if (!Array.isArray(skills) || !skills.length) return null;
    const byName = new Map(skills.map((skill) => [skill.name, skill]));
    return {
        name: READ_SKILL,
        keep: true,
        resultBudget: SKILL_RESULT_BUDGET,
        description: 'Read one of your skills — the method to follow for a kind of question (listed in the system '
            + 'prompt) — or, with `reference`, one of the reference files it lists. Read the skill before the analysis '
            + 'it covers, then follow its steps and its checklist.',
        parameters: {
            type: 'object',
            properties: {
                name: { type: 'string', enum: skills.map((skill) => skill.name), description: 'The skill.' },
                reference: { type: 'string', description: 'Optional: a reference the skill lists, by its name.' }
            },
            required: ['name']
        },
        label: ({ name, reference }) => (reference ? `Reading ${name}: ${reference}` : `Reading the ${name} skill`),
        run: async ({ name, reference }) => {
            const skill = byName.get(name);
            if (!skill) throw new Error(`no skill named "${name}": the skills are ${[...byName.keys()].join(', ')}`);
            const draftCheck = skill.draft
                ? [{ level: 'note', text: `${name} is a draft awaiting the user's interview: where it says "to confirm with the user", say which choice you made.` }]
                : [];
            if (reference != null && reference !== '') {
                const ref = skill.references.find((item) => item.name === reference);
                if (!ref) {
                    const names = skill.references.map((item) => item.name);
                    throw new Error(`${name} has no reference "${reference}"${names.length ? `; it lists ${names.join(', ')}` : ''}`);
                }
                if (typeof load !== 'function') throw new Error('the skill references are not available here');
                const text = await load(ref.path);
                return { data: { skill: name, reference: ref.name, text: String(text).trim() }, checks: draftCheck };
            }
            const refs = skill.references.length
                ? `\n\nReferences (read one with read_skill and \`reference\`):\n${skill.references.map((item) => `- ${item.name}: ${item.description}`).join('\n')}`
                : '';
            return {
                data: { skill: name, text: `${skill.body.replace(/^## References\s*\n[\s\S]*?(?=^## |(?![\s\S]))/m, '').trim()}${refs}` },
                checks: draftCheck
            };
        }
    };
};
