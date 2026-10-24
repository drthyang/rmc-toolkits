// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// The models the live arms (B, C) run, each through the same cases, prompt,
// tools and scoring, reported side by side. Local models cost nothing and
// may run dev and test; the API model runs the test set only, under a spend
// cap the tracker enforces (spend.js). Settings are the Copilot's own
// defaults (temperature 0.2, the turn's limits: runPolicy.js TURN_LIMITS),
// except where a model sets its own: `temperature` null sends no temperature
// at all, the warmer retry included, and `maxTokens` caps every reply
// (provider/model.js providerModel; the run policy records both). Only
// Claude Sonnet 5.5 sets them (the user's decisions, 2026-10-24); the local
// models' requests are as before.
//
// `server` is the server kind the model runs on, as the context guard names
// it (contextGuard.js detectServer: 'lmstudio', 'ollama', 'openai' for a
// cloud provider). The pre-run check (runPolicy.js preflight) refuses a run
// when the server there is another kind or does not answer, and an LM Studio
// model loaded with less than the run policy's minimum window (48k): the
// window is set when the model is loaded, not here.

import { MAX_REPLY_TOKENS } from './spend';

export const MODELS = Object.freeze({
    'qwen3.8-27b': {
        label: 'qwen3.8-27b (LM Studio)',
        baseUrl: 'http://localhost:1234/v1',
        model: 'qwen3.8-27b',
        server: 'lmstudio',
        api: false
    },
    'gemma4:26b': {
        label: 'gemma4:26b (Ollama)',
        baseUrl: 'http://localhost:11434/v1',
        model: 'gemma4:26b',
        server: 'ollama',
        api: false
    },
    'meta/muse-glimmer': {
        label: 'meta/muse-glimmer (LM Studio; 28B, tool use, up to 128k context)',
        baseUrl: 'http://localhost:1234/v1',
        model: 'meta/muse-glimmer',
        server: 'lmstudio',
        api: false
    },
    'claude-sonnet-5.5': {
        label: 'Claude Sonnet 5.5 (API)',
        baseUrl: 'https://api.anthropic.com/v1',
        model: 'claude-sonnet-5-5',
        server: 'openai',
        modelEnv: 'RMC_EVAL_SONNET_MODEL',
        apiKeyEnv: 'ANTHROPIC_API_KEY',
        api: true,
        // RMC's share of the cycle's $40 API total, which may be up to $28
        // (the user's decision, 2026-10-24: "Raise to $20").
        capUSD: 20,
        // Sonnet 5.5 takes only its default sampling: a request carrying a
        // temperature is refused (HTTP 400). None is sent, the warmer retry
        // included (null, as the sibling project NEBULA sets it; the user's
        // decision, 2026-10-24: "Fix both first").
        temperature: null,
        // Every reply capped at the reply the spend tracker prices
        // (max_tokens: hidden thinking, billed as output, counts within it),
        // so the tracker's worst case is a bound the API enforces (the same
        // decision).
        maxTokens: MAX_REPLY_TOKENS,
        splits: ['test'],
        runsPerCase: 3,
        // Arm B only: B and C together would roughly double the spend.
        arms: ['B'],
        // USD per million tokens (2026-10-14, from the orchestrator; checked
        // against the pricing page for claude-sonnet-5-5 on 2026-10-24). The
        // tracker charges every prompt token at the dearest prompt rate (a
        // cache write, $2.50), an upper bound whether or not a prompt is
        // cached: the client sees no usage report to tell.
        prices: { inputPerMTok: 2, outputPerMTok: 10, cacheReadPerMTok: 0.1, cacheWritePerMTok: 2.5 },
        priceEnv: { inputPerMTok: 'RMC_EVAL_PRICE_IN', outputPerMTok: 'RMC_EVAL_PRICE_OUT' }
    }
});

// The Copilot's temperature, for every model that does not set its own.
export const TEMPERATURE = 0.2;

/**
 * A model's settings with its environment overrides; throws when an API model
 * lacks its key, id or prices. `temperature` is the model's own when it sets
 * one (null: none sent), else TEMPERATURE; `maxTokens` only where set.
 */
export const resolveModel = (name, env = {}) => {
    const base = MODELS[name];
    if (!base) throw new Error(`unknown model '${name}'; known: ${Object.keys(MODELS).join(', ')}`);
    const model = { name, ...base, temperature: Object.hasOwn(base, 'temperature') ? base.temperature : TEMPERATURE };
    if (base.modelEnv && env[base.modelEnv]) model.model = env[base.modelEnv];
    if (base.priceEnv) {
        model.prices = {
            ...base.prices,
            inputPerMTok: Number(env[base.priceEnv.inputPerMTok] ?? base.prices.inputPerMTok),
            outputPerMTok: Number(env[base.priceEnv.outputPerMTok] ?? base.prices.outputPerMTok)
        };
    }
    if (base.apiKeyEnv) model.apiKey = env[base.apiKeyEnv] || null;
    if (model.api) {
        if (!model.model) throw new Error(`${name}: set ${base.modelEnv} to the exact model id`);
        if (!(model.prices.inputPerMTok > 0) || !(model.prices.outputPerMTok > 0)) {
            throw new Error(`${name}: set ${base.priceEnv.inputPerMTok} and ${base.priceEnv.outputPerMTok} (USD per million tokens) so the spend cap can be enforced`);
        }
        if (!model.apiKey) throw new Error(`${name}: set ${base.apiKeyEnv}`);
    }
    return model;
};
