import { useCallback, useEffect, useRef, useState } from 'react';
import { streamChat } from './provider/client';

// Shared streaming state for the summary, chat, and report views: accumulates
// deltas from streamChat into `text`, exposes stop(), and aborts on unmount.
// Stopping mid-stream keeps the partial text (it is still useful to read);
// only real failures surface as `error`.
export const useStreamedReply = () => {
    const [reply, setReply] = useState({ text: '', streaming: false, error: null });
    const abortRef = useRef(null);

    useEffect(() => () => abortRef.current?.abort(), []);

    const stop = useCallback(() => abortRef.current?.abort(), []);

    // Returns the full reply text, or null when the request failed.
    const start = useCallback(async ({ baseUrl, model, messages, temperature }) => {
        abortRef.current?.abort();
        const controller = new AbortController();
        abortRef.current = controller;
        setReply({ text: '', streaming: true, error: null });
        let text = '';
        try {
            for await (const delta of streamChat({ baseUrl, model, messages, temperature, signal: controller.signal })) {
                text += delta;
                setReply({ text, streaming: true, error: null });
            }
            setReply({ text, streaming: false, error: null });
            return text;
        } catch (error) {
            if (error?.name === 'AbortError') {
                setReply({ text, streaming: false, error: null });
                return text || null;
            }
            setReply({ text, streaming: false, error: error.message || 'The model request failed' });
            return null;
        }
    }, []);

    return { ...reply, start, stop };
};
