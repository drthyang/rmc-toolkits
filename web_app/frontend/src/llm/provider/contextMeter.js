// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// The context meter (THY-76): how full the active model's window is, read
// from the context guard's own numbers (contextGuard.js), so the meter and
// the guard never disagree. A request is the guard's estimate of its prompt
// (estimatePromptTokens) plus the reply reserve the guard keeps for that
// server; the window is the line the guard refuses at:
//   - Ollama: the model's maximum (POST /api/show), which the guard sizes
//     num_ctx up to and refuses above (estimate + REPLY_RESERVE_TOKENS > max);
//   - LM Studio: the context the model is loaded with (GET /api/v0/models),
//     refused below estimate + LMSTUDIO_REPLY_MIN_TOKENS;
//   - a cloud model: its published window (CLOUD_MODELS), with
//     REPLY_RESERVE_TOKENS of reply. The guard sends these unchanged; the
//     provider answers an over-long request with an error, never a cut.
// 100 % is exactly the guard's refusal line. The estimate errs high on
// purpose (2.5 characters a token, 1.4 for tool results and arguments;
// ~30 % above gemma4's real count on the Demo, README "The context window"),
// so the meter reads high, never low.

import {
    LMSTUDIO_REPLY_MIN_TOKENS, REPLY_RESERVE_TOKENS, detectServer, lmStudioLoadedContext, ollamaLoadedContext, ollamaMaxContext
} from './contextGuard';

// Neutral below 70 %, amber from 70 % to 90 %, red above 90 %.
export const AMBER_SHARE = 0.7;
export const RED_SHARE = 0.9;

// Published windows and prices (USD per million tokens) of the cloud models
// the meter knows, from Anthropic's model table (cached 2026-10-06). A model
// not listed has no window or price here: the meter shows its tokens only and
// says so. `tiers`: a price that changes with the prompt's size.
export const CLOUD_MODELS = Object.freeze({
    'claude-fable-5-1': { label: 'Claude Fable 5.1', window: 1_000_000, inputPerMTok: 10, outputPerMTok: 50 },
    'claude-fable-5': { label: 'Claude Fable 5', window: 1_000_000, inputPerMTok: 10, outputPerMTok: 50 },
    'claude-opus-5-5': { label: 'Claude Opus 5.5', window: 1_000_000, inputPerMTok: 4, outputPerMTok: 20 },
    'claude-opus-5': { label: 'Claude Opus 5', window: 1_000_000, inputPerMTok: 5, outputPerMTok: 25 },
    'claude-opus-4-8': { label: 'Claude Opus 4.8', window: 1_000_000, inputPerMTok: 5, outputPerMTok: 25 },
    'claude-opus-4-7': { label: 'Claude Opus 4.7', window: 1_000_000, inputPerMTok: 5, outputPerMTok: 25 },
    'claude-opus-4-6': { label: 'Claude Opus 4.6', window: 1_000_000, inputPerMTok: 5, outputPerMTok: 25 },
    'claude-sonnet-5-5': { label: 'Claude Sonnet 5.5', window: 1_000_000, inputPerMTok: 2, outputPerMTok: 10 },
    'claude-sonnet-5': { label: 'Claude Sonnet 5', window: 1_000_000, inputPerMTok: 2, outputPerMTok: 10 },
    'claude-sonnet-4-6': { label: 'Claude Sonnet 4.6', window: 1_000_000, inputPerMTok: 3, outputPerMTok: 15 },
    'claude-haiku-5-5': {
        label: 'Claude Haiku 5.5',
        window: 1_000_000,
        inputPerMTok: 0.1,
        outputPerMTok: 0.5,
        // Prompts above 100K tokens are priced higher.
        tiers: [{ above: 100_000, inputPerMTok: 0.5, outputPerMTok: 2.5 }]
    },
    'claude-haiku-4-5': { label: 'Claude Haiku 4.5', window: 200_000, inputPerMTok: 1, outputPerMTok: 5 }
});

/** A cloud model's entry by its id (a date suffix or "-latest" ignored), or null. */
export const cloudModel = (model) => {
    if (typeof model !== 'string' || !model) return null;
    const id = model.trim().toLowerCase().replace(/^.*\//, '').replace(/(@|-)(\d{8}|latest)$/, '');
    return CLOUD_MODELS[id] || null;
};

/** The meter's level for a share of the window (null: unknown window). */
export const meterLevel = (share) => {
    if (!Number.isFinite(share)) return 'neutral';
    if (share > RED_SHARE) return 'red';
    if (share >= AMBER_SHARE) return 'amber';
    return 'neutral';
};

/**
 * The window and reply reserve the guard decides a request against, from its
 * plan (planCopilotChat's, or probeLimits'): `{ server, window, reserve,
 * source, loadedContext }`. `window` null: not known (an Ollama model whose
 * maximum it does not report, LM Studio before the model is loaded, a cloud
 * model not in CLOUD_MODELS). `source` says where the window came from.
 */
export const meterLimits = (plan, { model = null } = {}) => {
    const server = plan?.server ?? 'openai';
    if (server === 'ollama') {
        const max = Number.isFinite(plan.maxContext) ? plan.maxContext : null;
        return {
            server, window: max, reserve: REPLY_RESERVE_TOKENS, source: max ? 'ollama_max' : null,
            loadedContext: Number.isFinite(plan.loadedContext) ? plan.loadedContext : null
        };
    }
    if (server === 'lmstudio') {
        const loaded = Number.isFinite(plan.loadedContext) ? plan.loadedContext : null;
        return { server, window: loaded, reserve: LMSTUDIO_REPLY_MIN_TOKENS, source: loaded ? 'lmstudio_loaded' : null, loadedContext: loaded };
    }
    const cloud = cloudModel(model);
    return { server, window: cloud?.window ?? null, reserve: REPLY_RESERVE_TOKENS, source: cloud ? 'published' : null, loadedContext: null };
};

/**
 * What the guard would plan against for this model, without a request: the
 * server kind and, for a local server, its windows (the same probes and
 * caches planCopilotChat uses; only localhost or Ollama's port is probed).
 */
export const probeLimits = async ({ baseUrl, model }, { signal } = {}) => {
    const detected = await detectServer(baseUrl, { signal });
    if (detected.kind === 'ollama') {
        const [loadedContext, maxContext] = await Promise.all([
            ollamaLoadedContext(baseUrl, model, { signal }),
            ollamaMaxContext(baseUrl, model, { signal })
        ]);
        return { server: 'ollama', loadedContext, maxContext };
    }
    if (detected.kind === 'lmstudio') {
        return { server: 'lmstudio', loadedContext: await lmStudioLoadedContext(baseUrl, model, { signal }) };
    }
    return { server: 'openai' };
};

/**
 * One request against its limits: `{ prompt, reserve, used, window, share,
 * level }`. used = prompt + reserve, the guard's `needed`; share = used /
 * window (null when the window is unknown); above 1 the guard refuses it.
 */
export const contextUsage = ({ estimate, limits }) => {
    const prompt = Math.max(0, Math.round(estimate || 0));
    const reserve = limits?.reserve ?? REPLY_RESERVE_TOKENS;
    const used = prompt + reserve;
    const window = Number.isFinite(limits?.window) && limits.window > 0 ? limits.window : null;
    const share = window ? used / window : null;
    return { prompt, reserve, used, window, share, level: meterLevel(share) };
};

/**
 * What a request costs on a priced cloud model (USD), from the same token
 * estimates: the prompt at the input price, the reply at the output price.
 * Null for a model with no price here (local models cost nothing).
 */
export const requestCost = ({ model, inputTokens, outputTokens = 0 }) => {
    const entry = cloudModel(model);
    if (!entry) return null;
    const tier = (entry.tiers || []).filter((item) => inputTokens > item.above).at(-1) || entry;
    return (inputTokens * tier.inputPerMTok + outputTokens * tier.outputPerMTok) / 1e6;
};
