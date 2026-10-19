// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// The run policy (plan item 1.2): the rules a live run ran under. The cases,
// the question, the prompt, the tools and the scoring are the same for every
// model; these are not all, so a comparison across models is partly a
// comparison of them. Each record keeps its policy (`record.policy`), each
// transcript prints it, and the report sets the columns' policies side by side.
//
// A policy states:
//   server         the kind the context guard detects (contextGuard.js
//                  detectServer: 'lmstudio', 'ollama', or 'openai' for every
//                  other server, the cloud providers among them), its version
//                  (Ollama), and the kind the model expects (models.js `server`);
//   temperature    the model's, and the floor of the request asked again
//                  after a looping reply (model.js ANTI_REPEAT_TEMPERATURE);
//   limits         the agent loop's budgets as the run used them
//                  (DEFAULT_LIMITS with the run's overrides);
//   context        the minimum window for this server kind (MIN_CONTEXT_TOKENS),
//                  how the window is sized, the window read from the server
//                  before the run (LM Studio's loaded context; Ollama's loaded
//                  context and the model's maximum) and after it (LM Studio),
//                  and every num_ctx an Ollama request was sent (from the
//                  guard's plan in each request as sent);
//   reply_cap      what bounds the reply: Ollama's num_predict; nothing is sent
//                  to LM Studio or a cloud provider;
//   context_guard  the guard's decisions, which the policy does not change: a
//                  request that would not fit is refused, never cut;
//   final_round    what the last request sends for tool_choice 'none', per
//                  server kind, as the model port sends it (pinned by the
//                  scripted check in __tests__/runPolicy.test.js);
//   history        one fresh question; earlier reasoning is never re-sent;
//   spend          an API model's cap and the reply each request is priced at;
//   warnings       what went wrong with the policy during the run.
//
// The decision of 2026-10-19 (the manager, on the user's delegation): the
// guard, its estimate (2.5 characters a token, 1.4 for tool JSON), the
// 8,192-token reply reserve, the budgets and the prompt stay as they are:
// refusing rather than cutting is the behaviour wanted. The live runs load
// each LM Studio model with a window of at least 48k, which holds the worst
// two-skill turn (src/__tests__/contextBudget.test.js), and the harness reads
// the loaded window before each run and refuses one below the minimum
// (`preflight`). Ollama sizes num_ctx itself (the guard) and is recorded only.

import { DEFAULT_LIMITS } from '../src/llm/copilot/agentLoop';
import {
    CHARS_PER_TOKEN, CONTEXT_STEP, JSON_CHARS_PER_TOKEN, LMSTUDIO_REPLY_MIN_TOKENS, MESSAGE_OVERHEAD_TOKENS, MIN_CONTEXT,
    READ_SHARE_LIMIT, REPLY_RESERVE_TOKENS, detectServer, lmStudioLoadedContext, ollamaLoadedContext, ollamaMaxContext
} from '../src/llm/provider/contextGuard';
import { ANTI_REPEAT_TEMPERATURE } from '../src/llm/provider/model';
import { MAX_REPLY_TOKENS } from './spend';

// The smallest window a live run may run under, per server kind (tokens).
// LM Studio fixes the window when the model is loaded: 48k (48 × 1,024 =
// 6 × CONTEXT_STEP) holds the worst two-skill turn of the dock and of arm B
// with the full reply reserve (contextBudget.test.js). Ollama: none, as the
// guard sizes num_ctx per request. A cloud provider's window is its own and
// is never probed.
export const MIN_CONTEXT_TOKENS = Object.freeze({ lmstudio: 49152, ollama: null, openai: null });

export const SERVER_LABELS = Object.freeze({
    lmstudio: 'LM Studio',
    ollama: 'Ollama',
    openai: 'an OpenAI-compatible server (cloud or other)'
});

// How the window is sized, per server kind.
export const CONTEXT_SIZING = Object.freeze({
    lmstudio: 'the window LM Studio loaded the model with, read before and after the run (GET /api/v0/models); '
        + 'the guard refuses a request it cannot hold with LMSTUDIO_REPLY_MIN_TOKENS of reply',
    ollama: 'num_ctx per request, sized by the guard: the loaded context when it holds the request and '
        + 'REPLY_RESERVE_TOKENS of reply, else the next multiple of CONTEXT_STEP (at least MIN_CONTEXT, '
        + 'at most the model\'s maximum); refused when even the maximum is too small',
    openai: 'the provider\'s own window; never probed, the request is sent unchanged'
});

// What bounds the reply, per server kind.
export const REPLY_CAPS = Object.freeze({
    lmstudio: Object.freeze({
        sent: null, rule: `none sent; the guard needs ${LMSTUDIO_REPLY_MIN_TOKENS.toLocaleString('en-US')} tokens of room beyond the prompt`
    }),
    ollama: Object.freeze({ sent: 'num_predict', rule: 'num_ctx minus the prompt estimate, per request' }),
    openai: Object.freeze({ sent: null, rule: 'none sent; the provider\'s own limit' })
});

// What the final round (the analysis budget used) sends for tool_choice
// 'none', per server kind. The model port decides the form
// (src/llm/provider/model.js toolChoiceForm, TOOL_CHOICE_FORMS): client.js
// puts `tool_choice` in every OpenAI-compatible body that has tools (OpenAI,
// Gemini, Anthropic, LM Studio), the tools still sent; Ollama's native
// /api/chat has no such field, so since protocol 2026-10-19.7 the port sends
// that request no tools, and the model cannot call one.
export const FINAL_ROUND_FORMS = Object.freeze({
    openai: Object.freeze({ tool_choice: 'none', tools: 'sent', how: '"tool_choice": "none" in the /chat/completions body (client.js)' }),
    lmstudio: Object.freeze({ tool_choice: 'none', tools: 'sent', how: '"tool_choice": "none" in the /v1/chat/completions body (client.js)' }),
    ollama: Object.freeze({
        tool_choice: null,
        tools: 'not sent',
        how: 'no tools: Ollama\'s native /api/chat has no tool_choice, so the model port enforces "none" by sending the '
            + 'final request without tools (provider/model.js); the earlier calls and results still travel as messages'
    })
});

// The context guard's decisions (provider/contextGuard.js), as constants.
export const CONTEXT_GUARD = Object.freeze({
    on_overflow: 'refuse',
    chars_per_token: CHARS_PER_TOKEN,
    json_chars_per_token: JSON_CHARS_PER_TOKEN,
    message_overhead_tokens: MESSAGE_OVERHEAD_TOKENS,
    reply_reserve_tokens: REPLY_RESERVE_TOKENS,
    lmstudio_reply_min_tokens: LMSTUDIO_REPLY_MIN_TOKENS,
    context_step: CONTEXT_STEP,
    min_context: MIN_CONTEXT,
    read_share_limit: READ_SHARE_LIMIT,
    probes: 'localhost, or port 11434, only'
});

/** The run policy refused a run before anything was sent. */
export class RunPolicyError extends Error {
    constructor(message, { server = null, needed = null } = {}) {
        super(message);
        this.name = 'RunPolicyError';
        this.server = server;
        this.needed = needed;
    }
}

const tokens = (count) => (Number.isFinite(count) ? Math.round(count).toLocaleString('en-US') : 'unknown');
const label = (kind) => SERVER_LABELS[kind] || String(kind);

/**
 * What the model's server is and the window it holds, read the way the
 * context guard reads it: `{ kind, version, loaded_tokens, max_tokens }`.
 * LM Studio: the loaded context (null when the model is not loaded under its
 * id). Ollama: the context the model is loaded with now and the model's
 * maximum. A cloud provider is not contacted (detectServer probes only
 * localhost and Ollama's port). Never throws, except on an abort.
 */
export const readServer = async (model, { signal } = {}) => {
    const detected = await detectServer(model.baseUrl, { signal });
    const server = { kind: detected.kind, version: detected.version ?? null, loaded_tokens: null, max_tokens: null };
    if (detected.kind === 'lmstudio') {
        server.loaded_tokens = await lmStudioLoadedContext(model.baseUrl, model.model, { signal });
    } else if (detected.kind === 'ollama') {
        const [loaded, max] = await Promise.all([
            ollamaLoadedContext(model.baseUrl, model.model, { signal }),
            ollamaMaxContext(model.baseUrl, model.model, { signal })
        ]);
        server.loaded_tokens = loaded;
        server.max_tokens = max;
    }
    return server;
};

const reloadHow = (id, needed) => `In LM Studio, eject ${id} and load it again with a Context Length of ${tokens(needed)} or more `
    + `(its load settings), or from a terminal: lms unload ${id} && lms load ${id} --context-length ${needed}. `
    + 'Then run the same command again: the runs already recorded are kept.';

/**
 * The pre-run check on what readServer found: the server is the kind the
 * model expects, and an LM Studio model is loaded with at least its minimum
 * window. Throws a RunPolicyError that says what was found, what is needed
 * and how to reload the model; returns `server` when the run may start.
 */
export const checkServer = (model, server, { minimum = MIN_CONTEXT_TOKENS } = {}) => {
    const expected = model.server ?? null;
    if (expected && server.kind !== expected) {
        throw new RunPolicyError(
            `${model.name} runs on ${label(expected)} at ${model.baseUrl}, but `
            + `${server.kind === 'openai' ? `no ${label(expected)} server answered there` : `the server there is ${label(server.kind)}`}, `
            + `so the run was not started. Start ${label(expected)}'s server${expected === 'lmstudio' ? ` and load ${model.model}` : ''}, `
            + 'then run the same command again.',
            { server }
        );
    }
    const needed = minimum[server.kind] ?? null;
    if (!Number.isFinite(needed)) return server;
    const id = model.model;
    if (!Number.isFinite(server.loaded_tokens)) {
        throw new RunPolicyError(
            `${label(server.kind)} has no model loaded as '${id}' (the id the harness sends), so its context cannot be read `
            + `and the run was not started: a request would load it at ${label(server.kind)}'s default context. The run policy `
            + `needs a context of at least ${tokens(needed)} tokens (48k). Load it: in LM Studio, load ${id} with a Context `
            + `Length of ${tokens(needed)} or more, or from a terminal: lms load ${id} --context-length ${needed} `
            + '(lms ps lists the loaded models and their ids). Then run the same command again: the runs already recorded are kept.',
            { server, needed }
        );
    }
    if (server.loaded_tokens < needed) {
        throw new RunPolicyError(
            `${label(server.kind)} runs ${id} with a ${tokens(server.loaded_tokens)}-token context, but the run policy needs at `
            + `least ${tokens(needed)} (48k) for an ${label(server.kind)} model, so the run was not started. ${reloadHow(id, needed)}`,
            { server, needed }
        );
    }
    return server;
};

/** The pre-run check: readServer, then checkServer. Resolves to the server it read; RunPolicyError below the policy. */
export const preflight = async (model, { minimum = MIN_CONTEXT_TOKENS, signal } = {}) => (
    checkServer(model, await readServer(model, { signal }), { minimum })
);

/**
 * A run's policy, as its record keeps it. `server` is what the pre-run check
 * read (null: not read, as in the scripted tests); `limits` the run's
 * overrides of the agent loop's budgets.
 */
export const buildRunPolicy = ({ model, limits = {}, server = null, minimum = MIN_CONTEXT_TOKENS }) => {
    const kind = server?.kind ?? null;
    return {
        server: { kind, expected: model?.server ?? null, version: server?.version ?? null, base_url: model?.baseUrl ?? null },
        temperature: { base: model?.temperature ?? null, warmer_retry_min: ANTI_REPEAT_TEMPERATURE },
        limits: { ...DEFAULT_LIMITS, ...limits },
        context: {
            min_tokens: kind ? (minimum[kind] ?? null) : null,
            min_tokens_by_server: { ...minimum },
            sizing: kind ? CONTEXT_SIZING[kind] : null,
            loaded_tokens: server?.loaded_tokens ?? null,
            model_max_tokens: server?.max_tokens ?? null,
            ...(kind === 'lmstudio' ? { loaded_tokens_after: null } : {}),
            num_ctx_used: []
        },
        reply_cap: kind ? { ...REPLY_CAPS[kind] } : null,
        context_guard: { ...CONTEXT_GUARD },
        final_round: kind ? { ...FINAL_ROUND_FORMS[kind] } : { tool_choice: 'none (as the loop asks)', tools: 'sent', how: 'server not read' },
        history: { earlier_turns: 0, reasoning_resent: false },
        spend: model?.api ? { cap_usd: model.capUSD ?? null, priced_reply_tokens: MAX_REPLY_TOKENS } : null,
        warnings: []
    };
};

/**
 * After the run: the num_ctx each request was sent joins its request
 * (`num_ctx`), read from the request as sent: the model port's prepare puts
 * the context guard's plan in it (`sent.plan.numCtx`, Ollama only; null on
 * other servers). A request the guard refused was never sent and has no
 * entry. The policy lists the values in order of first use.
 */
export const noteContextUsed = (record) => {
    if (!record?.policy) return record;
    record.requests.forEach((request) => {
        const numCtx = request.sent?.plan?.numCtx;
        request.num_ctx = Number.isFinite(numCtx) ? numCtx : null;
    });
    record.policy.context.num_ctx_used = [...new Set(record.requests.map((request) => request.num_ctx).filter(Number.isFinite))];
    return record;
};

/**
 * After the run: LM Studio's window read again. A window that changed (the
 * model reloaded, or unloaded and loaded again by a request at LM Studio's
 * default) is a warning: the run may not have run under the window read
 * before it.
 */
export const noteWindowAfter = (policy, after) => {
    if (!policy || policy.server.kind !== 'lmstudio') return policy;
    const now = after?.kind === 'lmstudio' && Number.isFinite(after.loaded_tokens) ? after.loaded_tokens : null;
    policy.context.loaded_tokens_after = now;
    if (now !== policy.context.loaded_tokens) {
        policy.warnings.push(`LM Studio's loaded context changed during the run: ${tokens(policy.context.loaded_tokens)} tokens before, `
            + `${now == null ? 'not loaded' : `${tokens(now)}`} after; the run may not have run under the window read before it`);
    }
    return policy;
};

// --- reading a policy ---------------------------------------------------------------

/** The window a run ran under, in a few words. */
export const windowSummary = (policy) => {
    if (!policy) return 'not recorded';
    const { server, context } = policy;
    if (!server.kind) return 'server not read';
    if (server.kind === 'lmstudio') {
        const after = 'loaded_tokens_after' in context && context.loaded_tokens_after != null && context.loaded_tokens_after !== context.loaded_tokens
            ? `, ${tokens(context.loaded_tokens_after)} after` : '';
        return `${tokens(context.loaded_tokens)} loaded${after} (minimum ${tokens(context.min_tokens)})`;
    }
    if (server.kind === 'ollama') {
        const used = context.num_ctx_used.length ? context.num_ctx_used.map(tokens).join(' → ') : 'none sent';
        return `num_ctx ${used} (sized per request; the model holds ${tokens(context.model_max_tokens)})`;
    }
    return 'the provider\'s own (not probed)';
};

export const limitsSummary = (limits) => `${limits.maxRounds} rounds × ${limits.maxCallsPerRound} calls, `
    + `${Math.round(limits.maxTurnMs / 1000)} s a turn, results ${tokens(limits.resultCharBudget)} / turn ${tokens(limits.turnCharBudget)} characters, `
    + `${limits.maxQuestions} questions, ${limits.maxEvidenceGates} gate send-backs`;

// The final round's form for a transcript line: tool_choice beside the tools
// (OpenAI-compatible servers), or no tools at all (Ollama since protocol
// 2026-10-19.7; an older Ollama record says "no tool_choice, tools still sent").
const finalRoundWords = (finalRound) => (finalRound.tools === 'not sent'
    ? 'no tools sent (the native API has no tool_choice)'
    : `${finalRound.tool_choice ? `tool_choice "${finalRound.tool_choice}"` : 'no tool_choice'}${finalRound.tools === 'sent' ? ', tools still sent' : ''}`);

/** The policy as transcript lines (Markdown list items). */
export const policyLines = (policy) => {
    if (!policy) return ['- Run policy: not recorded (a record from before protocol 2026-10-19.5)'];
    const { server, temperature, limits, reply_cap: replyCap, context_guard: guard, final_round: finalRound } = policy;
    const lines = [
        `- Server: ${server.kind ? label(server.kind) : 'not read (no pre-run check)'}${server.version ? ` ${server.version}` : ''}`
            + `${server.expected && server.expected !== server.kind ? ` (expected ${label(server.expected)})` : ''} · context window: ${windowSummary(policy)}`,
        `- Reply cap: ${replyCap ? replyCap.rule : 'server not read'} · final round: ${finalRoundWords(finalRound)}`,
        `- Temperature ${temperature.base} (a request asked again after a loop: ${temperature.warmer_retry_min} or more) · limits: ${limitsSummary(limits)}`,
        `- Context guard: ${guard.on_overflow === 'refuse' ? 'refuses a request that would not fit, never cuts it' : guard.on_overflow} `
            + `(${guard.chars_per_token} characters a token, ${guard.json_chars_per_token} for tool JSON; ${tokens(guard.reply_reserve_tokens)}-token reply reserve)`
            + `${policy.spend ? ` · spend cap $${policy.spend.cap_usd}` : ''}`
    ];
    (policy.warnings || []).forEach((warning) => lines.push(`- **Run policy warning:** ${warning}`));
    return lines;
};

/** The final round's form in a few words. */
export const finalRoundSummary = (policy) => {
    if (!policy) return 'not recorded';
    if (!policy.server.kind) return 'server not read';
    const { tool_choice: toolChoice, tools } = policy.final_round;
    // Ollama since protocol 2026-10-19.7: the port withholds the tools.
    if (tools === 'not sent') return 'no tools sent';
    return toolChoice ? `tool_choice "${toolChoice}" (tools sent)` : 'no tool_choice (tools sent)';
};
