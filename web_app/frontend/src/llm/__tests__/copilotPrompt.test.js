// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// The Copilot's prompt texts. How a turn assembles its messages from them —
// the context first, the prompt variant for the tools offered, the history
// window, the question last — is tested at the turn's interface
// (turn.test.js).

import { describe, expect, it } from 'vitest';
import {
    COPILOT_PROMPT,
    COPILOT_PROMPT_NO_TOOLS,
    COPILOT_PROMPT_VIEW,
    copilotSystemPrompt,
    toolBudgetNote
} from '../prompts/copilot';
import { SYSTEM_PROMPT } from '../prompts/system';

describe('the Copilot prompt', () => {
    it('extends the shared domain rules with the tool rules and the Outcome line', () => {
        expect(COPILOT_PROMPT.startsWith(SYSTEM_PROMPT)).toBe(true);
        expect(COPILOT_PROMPT).toContain('checks');
        // The model explains the Outcome; the app grades it (THY-18).
        expect(COPILOT_PROMPT).toContain('`Outcome: <one sentence>`');
        expect(COPILOT_PROMPT).toContain('Do not grade your answer');
        expect(COPILOT_PROMPT_NO_TOOLS).toContain('`Outcome: <one sentence>`');
        expect(COPILOT_PROMPT_NO_TOOLS).not.toContain('tools provided');
        // Numbers by reference and the evidence gate, with the tools only.
        expect(COPILOT_PROMPT).toContain('{{r2.peak.significance_sigma}}');
        expect(COPILOT_PROMPT).toContain('The app sends an answer back when its analysis did not run.');
        expect(COPILOT_PROMPT_NO_TOOLS).not.toContain('{{r2');
    });

    it('with the view tools, asks to point at cards and cite them as card: links', () => {
        expect(copilotSystemPrompt({ viewTools: true })).toBe(COPILOT_PROMPT_VIEW);
        expect(COPILOT_PROMPT_VIEW.startsWith(SYSTEM_PROMPT)).toBe(true);
        expect(COPILOT_PROMPT_VIEW).toContain('highlight_cards');
        expect(COPILOT_PROMPT_VIEW).toContain('[short label](card:<id>)');
        expect(COPILOT_PROMPT_VIEW).toContain('orient_3d_view');
        expect(COPILOT_PROMPT_VIEW).toContain('current_view');
        expect(COPILOT_PROMPT_VIEW).not.toContain('they never change what the user sees');
        expect(COPILOT_PROMPT).toContain('Tools only compute: they never change what the user sees.');
        // Tools off: no view rules either.
        expect(copilotSystemPrompt({ toolsEnabled: false, viewTools: true })).toBe(COPILOT_PROMPT_NO_TOOLS);
    });

    it('explains the budget note', () => {
        expect(toolBudgetNote(6)).toEqual({ role: 'user', content: expect.stringContaining('(6 rounds of tool calls)') });
    });
});
