import { useCallback, useEffect, useRef, useState } from 'react';
import { streamChat } from './provider/client';

// Streaming state for the chat view: accumulates the answer `text` and, for
// reasoning models, the separate `reasoning` stream plus `reasoningMs` (how
// long the model thought before the answer began). Exposes stop() and aborts on
// unmount. Stopping mid-stream keeps the partial output; only real failures
// surface as `error`.
const IDLE = { text: '', reasoning: '', reasoningMs: null, streaming: false, error: null };

export const useStreamedReply = () => {
    const [reply, setReply] = useState(IDLE);
    const abortRef = useRef(null);

    useEffect(() => () => abortRef.current?.abort(), []);

    const stop = useCallback(() => abortRef.current?.abort(), []);

    // Resolves to { text, reasoning, reasoningMs }, or null when the request
    // failed before producing anything.
    const start = useCallback(async ({ baseUrl, model, messages, temperature, apiKey }) => {
        abortRef.current?.abort();
        const controller = new AbortController();
        abortRef.current = controller;
        const startedAt = performance.now();
        setReply({ ...IDLE, streaming: true });
        let text = '';
        let reasoning = '';
        let reasoningMs = null;
        const settle = (extra) => {
            const result = { text, reasoning, reasoningMs, ...extra };
            setReply(result);
            return result;
        };
        try {
            for await (const chunk of streamChat({ baseUrl, model, messages, temperature, apiKey, signal: controller.signal })) {
                if (chunk.content) {
                    // First answer token: freeze how long the thinking phase took.
                    if (!text && reasoning) reasoningMs = performance.now() - startedAt;
                    text += chunk.content;
                }
                if (chunk.reasoning) reasoning += chunk.reasoning;
                setReply({ text, reasoning, reasoningMs, streaming: true, error: null });
            }
            const done = settle({ streaming: false, error: null });
            return { text: done.text, reasoning: done.reasoning, reasoningMs: done.reasoningMs };
        } catch (error) {
            if (error?.name === 'AbortError') {
                const done = settle({ streaming: false, error: null });
                return (done.text || done.reasoning) ? { text: done.text, reasoning: done.reasoning, reasoningMs: done.reasoningMs } : null;
            }
            settle({ streaming: false, error: error.message || 'The model request failed' });
            return null;
        }
    }, []);

    return { ...reply, start, stop };
};
