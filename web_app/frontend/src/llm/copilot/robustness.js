// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// What keeps a turn from hanging, looping or ending in nothing on a local
// model. The model the maintainer runs, qwen3.8-27b on LM Studio, is a dense
// reasoning model: 1–5 minutes a question in the 2026-10-10 spot check, once
// over 15 (one reasoning step took 10.5 minutes). gemma4:26b is faster but
// smaller. Each helper here is pure, or a thin wrapper, and tested alone:
//
//   - withDeadlines: a model reply that stalls (nothing for `stallMs`) or runs
//     past `maxMs` is cancelled and explained, never waited on forever;
//   - repeatsItself: a reply (or its reasoning) caught in a loop, writing the
//     same passage again and again, is stopped early instead of running until
//     the context is full;
//   - textToolCalls: a call a small model writes as JSON text instead of a
//     tool call is still run;
//   - keyNumbers: an older result compacted to fit the context keeps its
//     headline numbers, not only its checks;
//   - compactAnswer: a long earlier answer re-sent as memory keeps its
//     sentences with numbers and its verdict;
//   - explainModelError: an error the user can act on (out of memory, the
//     context overflowed, no model loaded, the server went away, a timeout).

export class ModelTimeoutError extends Error {
    constructor(kind, ms) {
        const minutes = Math.round(ms / 60000);
        const span = minutes >= 1 ? `${minutes} min` : `${Math.round(ms / 1000)} s`;
        super(kind === 'stall'
            ? `The model sent nothing for ${span}, so the request was stopped.`
            : `The model did not finish its reply within ${span}, so the request was stopped.`);
        this.name = 'ModelTimeoutError';
        this.kind = kind;
        this.ms = ms;
    }
}

/**
 * Iterate `stream` (an async iterator of model fragments) under two clocks:
 * no fragment for `stallMs`, or the whole reply past `maxMs`, calls `cancel()`
 * (which should abort the request) and throws a ModelTimeoutError. Reasoning
 * fragments count as progress for the stall clock, not for the total.
 */
export async function* withDeadlines(stream, { stallMs, maxMs, cancel = null }) {
    const iterator = stream[Symbol.asyncIterator]();
    const started = Date.now();
    let timer = null;
    let done = false;
    const deadline = (ms, kind) => new Promise((_, reject) => {
        timer = setTimeout(() => reject(new ModelTimeoutError(kind, kind === 'stall' ? stallMs : maxMs)), ms);
    });
    try {
        for (;;) {
            const left = Number.isFinite(maxMs) ? maxMs - (Date.now() - started) : Infinity;
            const stallLeft = Number.isFinite(stallMs) ? stallMs : Infinity;
            const limit = Math.min(left, stallLeft);
            let step;
            if (Number.isFinite(limit)) {
                const kind = left <= stallLeft ? 'max' : 'stall';
                try {
                    step = await Promise.race([iterator.next(), deadline(Math.max(0, limit), kind)]);
                } finally {
                    clearTimeout(timer);
                }
            } else {
                step = await iterator.next();
            }
            if (step.done) {
                done = true;
                return;
            }
            yield step.value;
        }
    } catch (error) {
        if (error instanceof ModelTimeoutError) cancel?.();
        throw error;
    } finally {
        // Stopped early (a deadline, or the reader broke off): close the stream.
        if (!done) Promise.resolve(iterator.return?.()).catch(() => {});
    }
}

/**
 * Whether `text` ends in a loop: its last `span` characters occur at least
 * `times` times in its last `reach` characters (a loop shorter than `span`
 * shows as overlapping copies). Real reasoning repeats short phrases, not
 * passages of this length.
 */
export const repeatsItself = (text, { span = 200, times = 4, reach = 8000 } = {}) => {
    const source = String(text || '');
    if (source.length < span) return false;
    const tail = source.slice(-span);
    if (!tail.trim()) return false;
    const region = source.slice(-reach);
    let count = 0;
    for (let at = region.indexOf(tail); at !== -1; at = region.indexOf(tail, at + 1)) {
        count += 1;
        if (count >= times) return true;
    }
    return false;
};

const fenced = (text) => {
    const match = /^```(?:json)?\s*([\s\S]*?)\s*```$/i.exec(text);
    return match ? match[1] : text;
};

/**
 * Tool calls written as text: the whole reply is one JSON object, or a list
 * of them, naming a known tool (`{"name": …, "arguments": {…}}`, or
 * `parameters` / `args`), optionally fenced or inside <tool_call> tags.
 * Returns [{ name, arguments }] (arguments as JSON text), or [] when the reply
 * is anything else — prose is never read as a call.
 */
export const textToolCalls = (text, toolNames) => {
    const known = new Set(toolNames || []);
    let source = String(text || '').trim();
    if (!source) return [];
    const tagged = [...source.matchAll(/<tool_call>\s*([\s\S]*?)\s*<\/tool_call>/g)];
    const pieces = tagged.length && !source.replace(/<tool_call>[\s\S]*?<\/tool_call>/g, '').trim()
        ? tagged.map((match) => match[1])
        : [fenced(source)];
    const calls = [];
    for (const piece of pieces) {
        let value;
        try {
            value = JSON.parse(piece);
        } catch {
            return [];
        }
        for (const item of Array.isArray(value) ? value : [value]) {
            const name = item?.name ?? item?.function?.name;
            const args = item?.arguments ?? item?.parameters ?? item?.args ?? item?.function?.arguments ?? {};
            if (typeof name !== 'string' || !known.has(name)) return [];
            if (typeof args === 'string') {
                calls.push({ name, arguments: args });
            } else if (args && typeof args === 'object' && !Array.isArray(args)) {
                calls.push({ name, arguments: JSON.stringify(args) });
            } else {
                return [];
            }
        }
    }
    return calls;
};

/**
 * The headline of a result for a compacted message: its top-level numbers,
 * short strings and booleans, and short lists of numbers (a window, a
 * direction), within `limit` characters of JSON. Larger fields are left out.
 */
export const keyNumbers = (data, limit = 600) => {
    if (!data || typeof data !== 'object' || Array.isArray(data)) return null;
    const out = {};
    let size = 2;
    for (const [key, value] of Object.entries(data)) {
        let kept;
        if (typeof value === 'number' && Number.isFinite(value)) kept = value;
        else if (typeof value === 'boolean') kept = value;
        else if (typeof value === 'string' && value.length <= 60 && key !== 'note') kept = value;
        else if (Array.isArray(value) && value.length <= 4 && value.every((item) => typeof item === 'number')) kept = value;
        if (kept === undefined) continue;
        const piece = JSON.stringify({ [key]: kept }).length - 1;
        if (size + piece > limit) break;
        out[key] = kept;
        size += piece;
    }
    return Object.keys(out).length ? out : null;
};

/**
 * An earlier answer shortened for the re-sent memory: under `limit`
 * characters it is unchanged; above, its sentences that hold a number, in
 * order, then the verdict line, up to the limit.
 */
export const compactAnswer = (text, limit = 2000) => {
    const answer = String(text || '');
    if (answer.length <= limit) return answer;
    const lines = answer.split('\n');
    const verdict = [...lines].reverse().find((line) => /^\W*outcome\b/i.test(line.trim())) || '';
    const body = verdict ? answer.slice(0, answer.lastIndexOf(verdict)) : answer;
    const sentences = body.replace(/\s+/g, ' ').split(/(?<=[.!?])\s+/).filter((sentence) => /\d/.test(sentence));
    const room = limit - verdict.length - 40;
    const kept = [];
    let size = 0;
    for (const sentence of sentences) {
        if (size + sentence.length + 1 > room) break;
        kept.push(sentence);
        size += sentence.length + 1;
    }
    return `${kept.join(' ')} (Shortened: the sentences with numbers are kept.)${verdict ? `\n${verdict.trim()}` : ''}`;
};

const PATTERNS = [
    [/out of memory|insufficient memory|OutOfMemory|failed to allocate|not enough memory/i,
        'The model server ran out of memory. Close other models (or other jobs using the GPU), load this one with a smaller context, or pick a smaller model.'],
    [/context (length|size|window)? ?(exceeded|overflow)|exceeds? (the )?(context|maximum context)|too many tokens|n_ctx|context overflow/i,
        'The conversation no longer fits the model\'s context. Start a new conversation, or load the model with a larger context.'],
    [/no models? (is |are )?loaded|model (is )?not loaded|model ['"]?[\w./:-]*['"]? not found|unknown model/i,
        'The model is not loaded on the server. Load it (LM Studio) or pull it (Ollama), then press Test in the settings.'],
    [/failed to fetch|networkerror|network error|load failed|ECONNREFUSED|connection (refused|reset|closed)/i,
        'The connection to the model server was lost. Check that it is still running.']
];

/**
 * A model error as the user should read it: the error's own message, and a
 * hint for what to do when it is one of the known failures. Returns
 * `{ message, hint }`.
 */
export const explainModelError = (error) => {
    const message = error?.message || String(error || 'The model request failed');
    if (error?.name === 'ModelTimeoutError') {
        return {
            message,
            hint: error.kind === 'stall'
                ? 'The server may be busy with another job or have stopped. Check it, then ask again.'
                : 'A reasoning model can think for minutes on one step. Ask a narrower question, or use a faster model (gemma4:26b).'
        };
    }
    const match = PATTERNS.find(([pattern]) => pattern.test(message));
    return { message, hint: match ? match[1] : null };
};
