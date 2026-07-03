// Minimal OpenAI-compatible chat client for local LLM servers (Ollama's /v1,
// LM Studio, or anything speaking the same dialect). Hand-rolled on fetch — no
// SDK, no SSE library — so the whole provider surface stays readable in one
// file and adds zero runtime dependencies.

const trimBase = (baseUrl) => (baseUrl || '').replace(/\/+$/, '');

// Cloud providers (OpenAI, Gemini, …) authenticate with a Bearer token; local
// servers need none, so the header is only added when a key is present.
const authHeaders = (apiKey) => (apiKey ? { Authorization: `Bearer ${apiKey}` } : {});

// A failed fetch to localhost surfaces as a bare TypeError both when the server
// is not running and when the browser blocked the response for CORS, so the
// hint has to name both causes — the user cannot tell them apart from the page.
const unreachableHint = (baseUrl) => (
    `Could not reach ${trimBase(baseUrl)}. Either the server is not running, or it is not `
    + 'allowing this site via CORS. For Ollama set OLLAMA_ORIGINS before `ollama serve`; '
    + 'for LM Studio enable CORS in the server settings. Some browsers (and future Chrome '
    + 'private-network rules) may also block public sites from calling localhost — if all '
    + 'else fails, run the app locally.'
);

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
    const response = await fetch(`${trimBase(baseUrl)}/models`, { signal, headers: authHeaders(apiKey) });
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
        headers: { 'Content-Type': 'application/json', ...authHeaders(apiKey) },
        body: JSON.stringify({ model, messages, temperature, stream }),
        signal
    });
    if (!response.ok) throw new Error(await describeHttpError(response));
    return response;
};

// Stream a chat completion, yielding content deltas as they arrive. The SSE
// body is `data: {json}` lines terminated by `data: [DONE]`; chunks can split
// mid-line, so incomplete tail lines are buffered across reads.
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
                const delta = parsed.choices?.[0]?.delta?.content;
                if (delta) yield delta;
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
