// Minimal OpenAI-compatible chat client for local LLM servers (Ollama's /v1,
// LM Studio, or anything speaking the same dialect). Hand-rolled on fetch — no
// SDK, no SSE library — so the whole provider surface stays readable in one
// file and adds zero runtime dependencies.

const trimBase = (baseUrl) => (baseUrl || '').replace(/\/+$/, '');

// Cloud providers (OpenAI, Gemini, …) authenticate with a Bearer token; local
// servers need none, so the header is only added when a key is present.
const authHeaders = (apiKey) => (apiKey ? { Authorization: `Bearer ${apiKey}` } : {});

// Anthropic's API blocks browser (CORS) requests unless the caller opts in with
// `anthropic-dangerous-direct-browser-access`, and its endpoints expect an API
// version. Its OpenAI-compatible /chat/completions still authenticates with a
// Bearer token, so only these two extra headers are Anthropic-specific.
const isAnthropic = (baseUrl) => {
    try {
        return new URL(baseUrl).hostname === 'api.anthropic.com';
    } catch {
        return false;
    }
};

// Headers for one request: Bearer auth plus, for Anthropic, its browser-access
// and version headers. Keyed off the base URL so no per-call plumbing is needed.
const requestHeaders = (baseUrl, apiKey) => ({
    ...authHeaders(apiKey),
    ...(isAnthropic(baseUrl)
        ? {
            'anthropic-version': '2023-06-01',
            'anthropic-dangerous-direct-browser-access': 'true'
        }
        : {})
});

// A failed fetch to localhost surfaces as a bare TypeError both when the server
// is not running and when the browser blocked the response for CORS, so the
// hint has to name both causes — the user cannot tell them apart from the page.
// It quotes this page's exact origin so the OLLAMA_ORIGINS value is copy-ready.
const unreachableHint = (baseUrl) => {
    const origin = typeof window !== 'undefined' ? window.location.origin : 'this page';
    return (
        `Could not reach ${trimBase(baseUrl)}. Either the server is not running, or it is not `
        + `allowing this page (${origin}) via CORS. Start Ollama with this origin allowed — `
        + `OLLAMA_ORIGINS="${origin}" ollama serve — or enable CORS in LM Studio's server settings. `
        + 'Safari also blocks HTTPS pages from calling http://localhost; use Chrome, Edge, or Firefox '
        + '(or run the app locally). See the setup guide for the macOS steps.'
    );
};

const describeHttpError = async (response) => {
    let detail = '';
    try {
        const payload = await response.json();
        detail = payload?.error?.message || payload?.error || '';
    } catch {
        // Non-JSON error bodies are fine; the status code is enough.
    }
    return `HTTP ${response.status}${detail ? `: ${detail}` : ''}`;
};

const httpHint = (status) => {
    if (status === 401) return 'Got 401 — the provider rejected the API key. Check that your key is correct and active.';
    if (status === 404) return 'Got 404 — check that the base URL ends in /v1 (e.g. http://localhost:11434/v1).';
    if (status === 403) return 'Got 403 — the server rejected this origin or key. Check its CORS/allowed-origins or key permissions.';
    if (status === 429) return 'Got 429 — the provider is rate-limiting or your quota is exhausted.';
    return null;
};

export const listModels = async (baseUrl, { signal, apiKey } = {}) => {
    const response = await fetch(`${trimBase(baseUrl)}/models`, { signal, headers: requestHeaders(baseUrl, apiKey) });
    if (!response.ok) {
        const error = new Error(await describeHttpError(response));
        error.status = response.status;
        throw error;
    }
    const payload = await response.json();
    return (payload.data || []).map((entry) => entry.id).filter(Boolean);
};

// Probe the server and translate failures into actionable setup hints.
// Returns { ok, models, error, hint } and never throws.
export const checkConnection = async (baseUrl, { signal, apiKey } = {}) => {
    try {
        const models = await listModels(baseUrl, { signal, apiKey });
        return { ok: true, models, error: null, hint: null };
    } catch (error) {
        if (error?.name === 'AbortError') throw error;
        const isNetworkError = error instanceof TypeError;
        return {
            ok: false,
            models: [],
            error: error.message || 'Connection failed',
            hint: isNetworkError ? unreachableHint(baseUrl) : httpHint(error.status)
        };
    }
};

const postChat = async ({ baseUrl, model, messages, temperature, stream, signal, apiKey }) => {
    const response = await fetch(`${trimBase(baseUrl)}/chat/completions`, {
        method: 'POST',
        headers: { 'Content-Type': 'application/json', ...requestHeaders(baseUrl, apiKey) },
        body: JSON.stringify({ model, messages, temperature, stream }),
        signal
    });
    if (!response.ok) throw new Error(await describeHttpError(response));
    return response;
};

// Stream a chat completion, yielding `{ content, reasoning }` deltas as they
// arrive (each chunk carries one or the other). The SSE body is `data: {json}`
// lines terminated by `data: [DONE]`; chunks can split mid-line, so incomplete
// tail lines are buffered across reads. Reasoning models (e.g. qwen3 via
// Ollama) stream their chain-of-thought in a separate `reasoning` field before
// the answer arrives in `content` — surfaced so the UI can show a thinking
// indicator instead of an empty bubble.
export async function* streamChat({ baseUrl, model, messages, temperature = 0.2, signal, apiKey }) {
    const response = await postChat({ baseUrl, model, messages, temperature, stream: true, signal, apiKey });
    if (!response.body) throw new Error('The server returned no response body to stream');
    const reader = response.body.getReader();
    const decoder = new TextDecoder();
    let buffer = '';
    try {
        for (;;) {
            const { done, value } = await reader.read();
            if (done) break;
            buffer += decoder.decode(value, { stream: true });
            const lines = buffer.split('\n');
            buffer = lines.pop() ?? '';
            for (const line of lines) {
                if (!line.startsWith('data:')) continue;
                const data = line.slice(5).trim();
                if (!data) continue;
                if (data === '[DONE]') return;
                let parsed;
                try {
                    parsed = JSON.parse(data);
                } catch {
                    continue;
                }
                const delta = parsed.choices?.[0]?.delta;
                if (!delta) continue;
                if (delta.content) yield { content: delta.content };
                // OpenAI-compatible servers name this `reasoning` (Ollama) or
                // `reasoning_content` (DeepSeek and others).
                const reasoning = delta.reasoning ?? delta.reasoning_content;
                if (reasoning) yield { reasoning };
            }
        }
    } finally {
        reader.cancel().catch(() => {});
    }
}

// Non-streaming completion, used where the whole reply is parsed at once
// (watchdog status lines, report narrative fallback).
export const completeChat = async ({ baseUrl, model, messages, temperature = 0, signal, apiKey }) => {
    const response = await postChat({ baseUrl, model, messages, temperature, stream: false, signal, apiKey });
    const payload = await response.json();
    const content = payload.choices?.[0]?.message?.content;
    if (typeof content !== 'string') throw new Error('The model returned no message content');
    return content;
};
