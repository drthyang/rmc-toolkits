// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// The Copilot's prompt texts. How a turn assembles its messages from them —
// the context first, the prompt variant for the tools offered, the history
// window, the question last — is tested at the turn's interface
// (turn.test.js).

import { describe, expect, it } from 'vitest';
import {
    ANSWER_RULES,
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

    // The user's report of 2026-11-01: answers overwhelmingly long, dumping
    // numbers. The answer's shape goes with the view tools (the dock, whose
    // step cards hold every number), before the Outcome line; the arms'
    // compute-only prompt (no view tools) is the one they always read.
    it('with the view tools, asks for the conclusion first and only the numbers that carry it', () => {
        expect(ANSWER_RULES).toMatch(/^Write the answer for the user reading the dock, and keep it short: first the conclusion/);
        expect(ANSWER_RULES).toContain('(at most five, unless the user asks for more)');
        expect(ANSWER_RULES).toContain('step card of each analysis');
        expect(COPILOT_PROMPT_VIEW).toContain(`${ANSWER_RULES} Finish every answer with one final line`);
        for (const method of ['full', 'off']) {
            expect(copilotSystemPrompt({ viewTools: true, askUser: true, method })).toContain(ANSWER_RULES);
            expect(copilotSystemPrompt({ viewTools: false, askUser: true, method })).not.toContain('reading the dock');
        }
        expect(COPILOT_PROMPT).not.toContain('reading the dock');
        expect(COPILOT_PROMPT_NO_TOOLS).not.toContain('reading the dock');
    });

    it('explains the budget note', () => {
        expect(toolBudgetNote(6)).toEqual({ role: 'user', content: expect.stringContaining('(6 rounds of tool calls)') });
    });
});
