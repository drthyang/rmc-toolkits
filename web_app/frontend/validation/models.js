// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// The models the live arms (B, C) run, each through the same cases, prompt,
// tools and scoring, reported side by side. Local models cost nothing and
// may run dev and test; the API model runs the test set only, under a spend
// cap the tracker enforces (spend.js). Settings are the Copilot's own
// defaults (temperature 0.2, the agent loop's DEFAULT_LIMITS).

export const MODELS = Object.freeze({
    'qwen3.8-27b': {
        label: 'qwen3.8-27b (LM Studio)',
        baseUrl: 'http://localhost:1234/v1',
        model: 'qwen3.8-27b',
        contextTokens: 32768,
        api: false
    },
    'gemma4:26b': {
        label: 'gemma4:26b (Ollama)',
        baseUrl: 'http://localhost:11434/v1',
        model: 'gemma4:26b',
        api: false
    },
    'meta/muse-glimmer': {
        label: 'meta/muse-glimmer (LM Studio; 28B, tool use, up to 128k context)',
        baseUrl: 'http://localhost:1234/v1',
        model: 'meta/muse-glimmer',
        contextTokens: 131072,
        api: false
    },
    'claude-sonnet-5.5': {
        label: 'Claude Sonnet 5.5 (API)',
        baseUrl: 'https://api.anthropic.com/v1',
        // The exact model id and the prices are set by whoever starts the run
        // (RMC_EVAL_SONNET_MODEL, RMC_EVAL_PRICE_IN, RMC_EVAL_PRICE_OUT in USD
        // per million tokens): the harness refuses an API run without them.
        model: null,
        modelEnv: 'RMC_EVAL_SONNET_MODEL',
        apiKeyEnv: 'ANTHROPIC_API_KEY',
        api: true,
        // RMC's share of the cycle's $40 API cap.
        capUSD: 12,
        splits: ['test'],
        runsPerCase: 3,
        // Arm B only by default: B and C together would roughly double the
        // spend (see the run-time and cost estimate in the validation doc).
        arms: ['B'],
        prices: { inputPerMTok: null, outputPerMTok: null },
        priceEnv: { inputPerMTok: 'RMC_EVAL_PRICE_IN', outputPerMTok: 'RMC_EVAL_PRICE_OUT' }
    }
});

export const TEMPERATURE = 0.2;

/** A model's settings with its environment overrides; throws when an API model lacks its id or prices. */
export const resolveModel = (name, env = {}) => {
    const base = MODELS[name];
    if (!base) throw new Error(`unknown model '${name}'; known: ${Object.keys(MODELS).join(', ')}`);
    const model = { name, ...base, temperature: TEMPERATURE };
    if (base.modelEnv && env[base.modelEnv]) model.model = env[base.modelEnv];
    if (base.priceEnv) {
        model.prices = {
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
