import { useEffect, useRef, useState } from 'react';
import { completeChat } from '../provider/client';
import { buildWatchdogMessages, parseWatchdogReply } from '../prompts/templates';
import { classifyConvergence, significantChange, watchdogStats } from './heuristics';

// Live convergence watchdog. No timers of its own: the existing Live Data poll
// already rebuilds rValueFile when the log files grow, so this hook just
// observes that prop. The heuristic classification is the source of truth and
// updates on every change for free; the LLM is only asked to narrate, and only
// when the status flips or enough genuinely new data arrived — never on the
// raw 3-second poll cadence.

const MS_PER_MINUTE = 60000;

export const useWatchdog = ({ rValueFile, settings }) => {
    const enabled = Boolean(settings?.watchdogEnabled);
    const history = rValueFile?.plotData?.series?.[0]?.y || null;
    // Key effects off content (length + last value), not object identity:
    // Live Data produces a fresh rValueFile object every poll even when the
    // underlying log did not change.
    const nSteps = history?.length || 0;
    const lastValue = nSteps ? history[nSteps - 1] : null;

    const [state, setState] = useState({
        status: 'unknown',
        source: 'heuristic',
        note: null,
        lastCheckedAt: null
    });
    // Latest history by ref (updated every render, before the analysis effect
    // below runs) so the effect can depend on cheap primitives instead of the
    // array identity, which Live Data churns every poll.
    const historyRef = useRef(null);
    useEffect(() => {
        historyRef.current = history;
    });
    const inFlightRef = useRef(false);
    const abortRef = useRef(null);
    // Last LLM call bookkeeping for the throttle gate.
    const lastLlmRef = useRef({ status: null, stats: null, at: 0 });

    useEffect(() => () => abortRef.current?.abort(), []);

    const { baseUrl, model, watchdogIntervalMin = 5 } = settings || {};

    useEffect(() => {
        if (!enabled || nSteps < 2) return;
        const values = historyRef.current;
        const status = classifyConvergence(values);
        const stats = watchdogStats(values);

        setState((current) => (
            current.status === status
                ? current
                : { ...current, status, source: 'heuristic', note: null }
        ));

        if (!baseUrl || !model || inFlightRef.current) return;
        const statusChanged = status !== lastLlmRef.current.status;
        const intervalElapsed = Date.now() - lastLlmRef.current.at >= watchdogIntervalMin * MS_PER_MINUTE;
        if (!statusChanged && !(intervalElapsed && significantChange(lastLlmRef.current.stats, stats))) {
            return;
        }

        inFlightRef.current = true;
        const controller = new AbortController();
        abortRef.current = controller;
        completeChat({
            baseUrl,
            model,
            messages: buildWatchdogMessages(stats, status, lastLlmRef.current.status),
            temperature: 0,
            signal: controller.signal
        }).then((reply) => {
            const parsed = parseWatchdogReply(reply);
            setState((current) => ({
                ...current,
                source: 'llm',
                note: parsed?.note || reply.trim().slice(0, 200) || null,
                lastCheckedAt: Date.now()
            }));
        }).catch(() => {
            // Unreachable server or a malformed reply: the heuristic badge is
            // already correct, so failures are silent by design.
        }).finally(() => {
            // Record the attempt either way so a broken server is retried on
            // the throttle schedule instead of every poll.
            lastLlmRef.current = { status, stats, at: Date.now() };
            inFlightRef.current = false;
        });
    }, [enabled, nSteps, lastValue, baseUrl, model, watchdogIntervalMin]);

    if (!enabled || nSteps < 2) {
        return { status: 'off', source: 'heuristic', note: null, lastCheckedAt: null };
    }
    return state;
};
