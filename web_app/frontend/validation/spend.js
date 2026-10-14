// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// The API spend tracker. Every request to a paid model is priced before it
// is sent (its prompt estimate plus the largest reply the loop allows) and
// refused when that could cross the cap; after it, the estimate of what was
// used is charged to a ledger that persists across processes, so a restarted
// batch cannot overspend. The client does not report token usage, so the
// charge is the Copilot's own conservative estimate (contextGuard: 2.5
// characters per token of text, 1.4 of JSON), an over-estimate for Claude.

import { CHARS_PER_TOKEN, estimatePromptTokens } from '../src/llm/provider/contextGuard';

export class SpendCapError extends Error {
    constructor(message) {
        super(message);
        this.name = 'SpendCapError';
    }
}

// The reply budget a single request may use (tokens); the price check
// reserves it before sending.
export const MAX_REPLY_TOKENS = 8192;

export const createSpendTracker = ({ capUSD, prices, ledger }) => {
    if (!(capUSD > 0)) throw new Error('a spend cap in USD is required');
    if (!(prices?.inputPerMTok > 0) || !(prices?.outputPerMTok > 0)) throw new Error('prices per million tokens are required');
    // Prompt tokens at the dearest prompt rate: a cache write can cost more
    // than plain input, and the client cannot tell which was billed.
    const promptPerMTok = Math.max(prices.inputPerMTok, prices.cacheWritePerMTok || 0, prices.cacheReadPerMTok || 0);
    const cost = (inputTokens, outputTokens) => (inputTokens * promptPerMTok + outputTokens * prices.outputPerMTok) / 1e6;
    const spent = () => ledger.load().reduce((sum, entry) => sum + (entry.usd || 0), 0);
    return {
        capUSD,
        spent,
        cost,
        /** Throws SpendCapError when this request could cross the cap. */
        reserve: (request) => {
            const inputTokens = estimatePromptTokens({ messages: request.messages, tools: request.tools ?? null });
            const worst = cost(inputTokens, MAX_REPLY_TOKENS);
            const now = spent();
            if (now + worst > capUSD) {
                throw new SpendCapError(`spend cap reached: $${now.toFixed(2)} spent of $${capUSD.toFixed(2)}; `
                    + `the next request could cost up to $${worst.toFixed(2)}`);
            }
            return { inputTokens };
        },
        /** Records what a request used (estimated from what was streamed). */
        charge: ({ inputTokens, outputChars, note }) => {
            const outputTokens = Math.ceil(outputChars / CHARS_PER_TOKEN);
            const usd = cost(inputTokens, outputTokens);
            ledger.append({ at: new Date().toISOString(), note, inputTokens, outputTokens, usd: Number(usd.toFixed(6)) });
            return usd;
        }
    };
};

/** An in-memory ledger (tests); the CLI uses a JSON-lines file. */
export const memoryLedger = (entries = []) => ({ load: () => entries, append: (entry) => { entries.push(entry); } });

/**
 * A model call (contextGuard's streamCopilotChat signature) that reserves
 * before each request and charges after it, counting what streamed back.
 */
export const meteredCall = (callModel, tracker, note) => async function* metered(request) {
    const { inputTokens } = tracker.reserve(request);
    let outputChars = 0;
    try {
        for await (const chunk of callModel(request)) {
            outputChars += (chunk.content?.length || 0) + (chunk.reasoning?.length || 0)
                + (chunk.toolCall?.arguments?.length || 0) + (chunk.toolCall?.name?.length || 0);
            yield chunk;
        }
    } finally {
        tracker.charge({ inputTokens, outputChars, note });
    }
};
