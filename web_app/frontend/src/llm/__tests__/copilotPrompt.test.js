// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

import { describe, expect, it } from 'vitest';
import {
    CHAT_HISTORY_TURNS,
    COPILOT_PROMPT,
    COPILOT_PROMPT_NO_TOOLS,
    COPILOT_PROMPT_VIEW,
    buildCopilotMessages,
    toolBudgetNote
} from '../prompts/copilot';
import { SYSTEM_PROMPT } from '../prompts/system';
import { conversationMemory } from '../useCopilot';

const context = { run: 'test-run', live_mode: false };

describe('buildCopilotMessages', () => {
    it('injects the context before the conversation and appends the question', () => {
        const messages = buildCopilotMessages({ context, history: [], userText: 'which dataset fits worst?' });
        expect(messages[0]).toEqual({ role: 'system', content: COPILOT_PROMPT });
        expect(messages[1].content).toContain('"run": "test-run"');
        expect(messages.at(-1)).toEqual({ role: 'user', content: 'which dataset fits worst?' });
    });

    it('keeps only the newest turns', () => {
        const history = Array.from({ length: 20 }, (_, index) => ({
            role: index % 2 ? 'assistant' : 'user',
            content: `turn ${index}`
        }));
        const messages = buildCopilotMessages({ context, history, userText: 'q' });
        const historyMessages = messages.slice(3, -1);
        expect(historyMessages).toHaveLength(CHAT_HISTORY_TURNS);
        expect(historyMessages[0].content).toBe(`turn ${20 - CHAT_HISTORY_TURNS}`);
    });

    it('uses the tool-free prompt when tools are off', () => {
        const messages = buildCopilotMessages({ context, userText: 'q', toolsEnabled: false });
        expect(messages[0].content).toBe(COPILOT_PROMPT_NO_TOOLS);
        expect(messages[2].content).not.toContain('analyses I run');
    });
});

describe('the Copilot prompt', () => {
    it('extends the shared domain rules with the tool rules and the Outcome line', () => {
        expect(COPILOT_PROMPT.startsWith(SYSTEM_PROMPT)).toBe(true);
        expect(COPILOT_PROMPT).toContain('checks');
        expect(COPILOT_PROMPT).toContain('Outcome: achieved|partial|not achieved');
        expect(COPILOT_PROMPT_NO_TOOLS).toContain('Outcome: achieved|partial|not achieved');
        expect(COPILOT_PROMPT_NO_TOOLS).not.toContain('tools provided');
    });

    it('with the view tools, asks to point at cards and cite them as card: links', () => {
        const messages = buildCopilotMessages({ context, userText: 'q', viewTools: true });
        expect(messages[0].content).toBe(COPILOT_PROMPT_VIEW);
        expect(COPILOT_PROMPT_VIEW.startsWith(SYSTEM_PROMPT)).toBe(true);
        expect(COPILOT_PROMPT_VIEW).toContain('highlight_cards');
        expect(COPILOT_PROMPT_VIEW).toContain('[short label](card:<id>)');
        expect(COPILOT_PROMPT_VIEW).toContain('orient_3d_view');
        expect(COPILOT_PROMPT_VIEW).toContain('current_view');
        expect(COPILOT_PROMPT_VIEW).not.toContain('they never change what the user sees');
        expect(COPILOT_PROMPT).toContain('Tools only compute: they never change what the user sees.');
        // Tools off: no view rules either.
        expect(buildCopilotMessages({ context, userText: 'q', toolsEnabled: false, viewTools: true })[0].content)
            .toBe(COPILOT_PROMPT_NO_TOOLS);
    });

    it('explains the budget note', () => {
        expect(toolBudgetNote(6)).toEqual({ role: 'user', content: expect.stringContaining('(6 rounds of tool calls)') });
    });
});

describe('conversationMemory', () => {
    it('pairs each answered question with its answer and drops unanswered ones', () => {
        const turns = [
            { role: 'user', content: 'q1' },
            { role: 'assistant', text: 'a1\nOutcome: achieved — x' },
            { role: 'user', content: 'q2 (stopped before an answer)' },
            { role: 'assistant', text: '   ' },
            { role: 'user', content: 'q3 (failed: no assistant turn)' },
            { role: 'user', content: 'q4' },
            { role: 'assistant', text: 'a4' }
        ];
        expect(conversationMemory(turns)).toEqual([
            { role: 'user', content: 'q1' },
            { role: 'assistant', content: 'a1\nOutcome: achieved — x' },
            { role: 'user', content: 'q4' },
            { role: 'assistant', content: 'a4' }
        ]);
    });
});
