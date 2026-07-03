import { describe, expect, it } from 'vitest';
import {
    CHAT_HISTORY_TURNS,
    buildChatMessages,
    buildWatchdogMessages,
    parseWatchdogReply
} from '../prompts/templates';

const context = { run: 'test-run', live_mode: false };

describe('buildChatMessages', () => {
    it('injects context before the conversation and appends the question', () => {
        const messages = buildChatMessages(context, [], 'which dataset fits worst?');
        expect(messages[0].role).toBe('system');
        expect(messages[1].content).toContain('"run": "test-run"');
        expect(messages[messages.length - 1]).toEqual({ role: 'user', content: 'which dataset fits worst?' });
    });

    it('keeps only the newest turns', () => {
        const history = Array.from({ length: 20 }, (_, index) => ({
            role: index % 2 ? 'assistant' : 'user',
            content: `turn ${index}`
        }));
        const messages = buildChatMessages(context, history, 'q');
        const historyMessages = messages.slice(3, -1);
        expect(historyMessages).toHaveLength(CHAT_HISTORY_TURNS);
        expect(historyMessages[0].content).toBe(`turn ${20 - CHAT_HISTORY_TURNS}`);
    });
});

describe('watchdog prompt round trip', () => {
    it('mentions the stats and the heuristic label', () => {
        const stats = { n_steps: 500, last: 0.87 };
        const messages = buildWatchdogMessages(stats, 'improving', 'stalled');
        expect(messages[1].content).toContain('"n_steps":500');
        expect(messages[1].content).toContain('"improving"');
        expect(messages[1].content).toContain('"stalled"');
    });

    it('parses well-formed replies', () => {
        expect(parseWatchdogReply('STATUS: improving — chi dropped to 0.87')).toEqual({
            status: 'improving',
            note: 'chi dropped to 0.87'
        });
        expect(parseWatchdogReply('status: DIVERGING - rising since step 400').status).toBe('diverging');
        expect(parseWatchdogReply('STATUS: converged')).toEqual({ status: 'converged', note: null });
    });

    it('returns null on replies that ignore the format', () => {
        expect(parseWatchdogReply('The refinement looks fine to me.')).toBeNull();
        expect(parseWatchdogReply('')).toBeNull();
        expect(parseWatchdogReply(null)).toBeNull();
    });
});
