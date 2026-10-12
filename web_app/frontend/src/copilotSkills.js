// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// The AI Copilot's skills, the host half of the llm module's `skills` prop:
// the repository's Claude Code skills (`.claude/skills/<name>/SKILL.md` at the
// repo root), so the in-app agent and Claude Code read the same methods. The
// SKILL.md files ship with the app; the reference files beside them are
// loaded only when the agent reads one (read_skill with `reference`). The
// module parses them (llm/copilot/skills.js).

// Paths as the repository names them: '../../../.claude/…' → '.claude/…'.
const repoPath = (key) => key.replace(/^(\.\.\/)+/, '');

const SKILL_FILES = import.meta.glob('../../../.claude/skills/*/SKILL.md', { query: '?raw', import: 'default', eager: true });
const REFERENCE_FILES = import.meta.glob(
    ['../../../.claude/skills/*/*.md', '!../../../.claude/skills/*/SKILL.md'],
    { query: '?raw', import: 'default' }
);
const referenceLoaders = new Map(Object.entries(REFERENCE_FILES).map(([key, load]) => [repoPath(key), load]));

export const copilotSkills = {
    files: Object.fromEntries(Object.entries(SKILL_FILES).map(([key, text]) => [repoPath(key), text])),
    references: [...referenceLoaders.keys()],
    load: async (path) => {
        const load = referenceLoaders.get(path);
        if (!load) throw new Error(`no bundled skill file ${path}`);
        return load();
    }
};
