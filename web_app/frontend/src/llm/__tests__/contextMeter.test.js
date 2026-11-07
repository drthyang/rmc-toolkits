// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// The context meter reads the guard's own numbers (THY-76): 100 % is exactly
// where the guard refuses, on Ollama and on LM Studio; the levels' bounds;
// the windows and prices of the cloud models it knows; and the request it
// previews is the one the turn sends first.

import { afterEach, beforeEach, describe, expect, it, vi } from 'vitest';
import {
    ContextLimitError, LMSTUDIO_REPLY_MIN_TOKENS, REPLY_RESERVE_TOKENS, estimatePromptTokens, forgetServers, planContext, planCopilotChat
} from '../provider/contextGuard';
import {
    AMBER_SHARE, RED_SHARE, cloudModel, contextUsage, meterLevel, meterLimits, probeLimits, requestCost
} from '../provider/contextMeter';
import { historyWindow, requestPreview, runTurn } from '../copilot/turn';
import { HISTORY_MESSAGES } from '../prompts/copilot';

const json = (body, status = 200) => ({ ok: status < 400, status, json: async () => body });

describe('the meter\'s levels', () => {
    it('is neutral below 70 %, amber from 70 % to 90 %, red above', () => {
        expect([AMBER_SHARE, RED_SHARE]).toEqual([0.7, 0.9]);
        expect(meterLevel(0.69)).toBe('neutral');
        expect(meterLevel(0.7)).toBe('amber');
        expect(meterLevel(0.9)).toBe('amber');
        expect(meterLevel(0.901)).toBe('red');
        expect(meterLevel(1.4)).toBe('red');
        // An unknown window has no share: no warning colour from a guess.
        expect(meterLevel(null)).toBe('neutral');
    });

    it('counts the reply reserve: used = the prompt + the reserve, against the window', () => {
        const usage = contextUsage({ estimate: 12000, limits: { window: 32768, reserve: REPLY_RESERVE_TOKENS } });
        expect(usage).toEqual({ prompt: 12000, reserve: 8192, used: 20192, window: 32768, share: 20192 / 32768, level: 'neutral' });
        expect(contextUsage({ estimate: 12000, limits: { window: null, reserve: 8192 } })).toMatchObject({ used: 20192, share: null, level: 'neutral' });
    });
});

describe('100 % is the guard\'s refusal line', () => {
    it('Ollama: the model\'s maximum, with the guard\'s reply reserve', () => {
        const limits = meterLimits({ server: 'ollama', maxContext: 40960, loadedContext: 16384 });
        expect(limits).toMatchObject({ server: 'ollama', window: 40960, reserve: REPLY_RESERVE_TOKENS, source: 'ollama_max', loadedContext: 16384 });
        for (let estimate = 30000; estimate <= 34000; estimate += 97) {
            const refused = Boolean(planContext({ estimate, max: 40960 }).refuse);
            expect(contextUsage({ estimate, limits }).share > 1).toBe(refused);
        }
        // No maximum reported: the guard never refuses, the meter knows no window.
        expect(meterLimits({ server: 'ollama', maxContext: null })).toMatchObject({ window: null, source: null });
    });

    describe('LM Studio: the loaded context, with its 4,096 tokens of reply', () => {
        let loaded;
        beforeEach(() => {
            forgetServers();
            vi.stubGlobal('fetch', vi.fn(async (url) => {
                const path = new URL(url).pathname;
                if (path === '/api/v0/models') return json({ data: [{ id: 'm', state: 'loaded', loaded_context_length: loaded }] });
                return json({ error: 'not found' }, 404);
            }));
        });
        afterEach(() => vi.unstubAllGlobals());

        it('refuses exactly where the meter passes 100 %', async () => {
            const request = { baseUrl: 'http://localhost:1234/v1', model: 'm', messages: [{ role: 'user', content: 'x'.repeat(50000) }] };
            const estimate = estimatePromptTokens(request);
            for (const window of [estimate + LMSTUDIO_REPLY_MIN_TOKENS - 1, estimate + LMSTUDIO_REPLY_MIN_TOKENS, estimate + 9000]) {
                loaded = window;
                const limits = meterLimits(await probeLimits(request), { model: 'm' });
                expect(limits).toMatchObject({ server: 'lmstudio', window, reserve: LMSTUDIO_REPLY_MIN_TOKENS, source: 'lmstudio_loaded' });
                const over = contextUsage({ estimate, limits }).share > 1;
                const outcome = await planCopilotChat(request).then(() => 'sent', (error) => error);
                if (over) expect(outcome).toBeInstanceOf(ContextLimitError);
                else expect(outcome).toBe('sent');
            }
        });
    });
});

describe('cloud models', () => {
    it('knows the published windows and prices of the Claude models, by id', () => {
        expect(cloudModel('claude-sonnet-5-5')).toMatchObject({ window: 1_000_000, inputPerMTok: 2, outputPerMTok: 10 });
        expect(cloudModel('claude-haiku-4-5-20251001')).toMatchObject({ window: 200_000 });
        expect(cloudModel('Claude-Opus-5-5')).toMatchObject({ inputPerMTok: 4, outputPerMTok: 20 });
        expect(cloudModel('gpt-unknown')).toBeNull();
        expect(meterLimits({ server: 'openai' }, { model: 'claude-sonnet-5-5' }))
            .toMatchObject({ window: 1_000_000, reserve: REPLY_RESERVE_TOKENS, source: 'published' });
        // A model it does not know: the tokens only, no window.
        expect(meterLimits({ server: 'openai' }, { model: 'gpt-x' })).toMatchObject({ window: null, source: null });
    });

    it('prices a request from the same token counts; local and unknown models have no price', () => {
        expect(requestCost({ model: 'claude-sonnet-5-5', inputTokens: 1_000_000 })).toBeCloseTo(2);
        expect(requestCost({ model: 'claude-sonnet-5-5', inputTokens: 10000, outputTokens: 1000 })).toBeCloseTo(0.03);
        // Haiku 5.5: prompts above 100K tokens at the higher price.
        expect(requestCost({ model: 'claude-haiku-5-5', inputTokens: 100_000 })).toBeCloseTo(0.01);
        expect(requestCost({ model: 'claude-haiku-5-5', inputTokens: 200_000 })).toBeCloseTo(0.1);
        expect(requestCost({ model: 'gemma4:26b', inputTokens: 1000 })).toBeNull();
    });

    it('never probes a cloud server for its window', async () => {
        const fetchMock = vi.fn();
        vi.stubGlobal('fetch', fetchMock);
        await expect(probeLimits({ baseUrl: 'https://api.anthropic.com/v1', model: 'claude-sonnet-5-5' })).resolves.toEqual({ server: 'openai' });
        expect(fetchMock).not.toHaveBeenCalled();
        vi.unstubAllGlobals();
    });
});

describe('the request the meter previews', () => {
    const run = {
        runName: 'demo', plotFiles: [], structure: null, symmetry: null, runSettings: null, chiLog: null, watchdog: null, extras: null
    };
    const tool = {
        name: 'compute_bond_angles',
        description: 'Bond angles.',
        parameters: { type: 'object', properties: { apex: { type: 'string' } } },
        run: async () => ({ data: { n: 1 }, checks: [] })
    };
    const history = [
        { question: 'Has it converged?', answer: 'Not yet: χ² still falls 2 % per 1000 moves.\nOutcome: achieved — falling.' },
        { question: 'Stopped one', answer: '' }
    ];

    it('is the first request the turn sends: the same messages and tools', async () => {
        const spec = { question: 'What are the angles?', run, tools: [tool], askUser: async () => ({ source: 'skipped' }), history };
        const preview = requestPreview(spec);
        let sent = null;
        for await (const event of runTurn({
            ...spec,
            model: async function* () { yield { content: 'Done.\nOutcome: achieved — x.' }; }
        })) {
            if (event.type === 'request' && !sent) sent = event.sent;
        }
        expect(preview.messages).toEqual(sent.messages);
        expect(preview.tools).toEqual(sent.tools);
        expect(estimatePromptTokens(preview)).toBe(estimatePromptTokens(sent));
        // Tools off: no tools in the request.
        expect(requestPreview({ ...spec, config: { toolsEnabled: false } }).tools).toBeNull();
    });

    it('says what the history window leaves out: the oldest questions, the long answers', () => {
        const keep = HISTORY_MESSAGES / 2;
        const many = Array.from({ length: keep + 2 }, (_, index) => ({
            question: `q${index}`, answer: index === keep + 1 ? `${'Long. '.repeat(400)}` : `a${index}`
        }));
        expect(historyWindow(many)).toEqual({ answered: keep + 2, sent: keep, dropped: 2, shortened: 1 });
        expect(historyWindow([...history])).toEqual({ answered: 1, sent: 1, dropped: 0, shortened: 0 });
        expect(historyWindow([])).toEqual({ answered: 0, sent: 0, dropped: 0, shortened: 0 });
    });
});
