// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

import { useEffect, useRef, useState } from 'react';
import { completeChat } from '../provider/client';
import { temperatureFor } from '../provider/sampling';
import { buildWatchdogMessages, parseWatchdogReply } from '../prompts/templates';
import { classifyRun, significantChange, watchdogRunStats } from './heuristics';

// Live convergence watchdog. No timers of its own: the existing Live Data poll
// already rebuilds rValueFile when the log files grow, so this hook just
// observes that prop. The heuristic classification is the source of truth and
// updates on every change for free: every chi^2 term of the .log and their
// total, by the rule the Copilot's assess_convergence reads (heuristics.js
// judgeConvergence; the chart draws the last column only). The LLM is only
// asked to narrate, and only when the status flips or enough genuinely new
// data arrived — never on the raw 3-second poll cadence. Its request runs at
// temperature 0, or with none for a model that takes none (Claude:
// provider/sampling.js).

const MS_PER_MINUTE = 60000;

export const useWatchdog = ({ rValueFile, settings }) => {
    const enabled = Boolean(settings?.watchdogEnabled);
    const plotData = rValueFile?.plotData || null;
    const history = plotData?.series?.[0]?.y || null;
    // The .log column the chart draws (one fit term's chi^2, not the total),
    // and whether the verdict judged every term and their total (a log that
    // names its chi^2 columns) or that column alone.
    const column = plotData?.chiColumn || null;
    const allTerms = Array.isArray(plotData?.chiTerms?.total);
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
    // Latest plot data by ref (updated every render, before the analysis
    // effect below runs) so the effect can depend on cheap primitives instead
    // of the array identity, which Live Data churns every poll.
    const plotDataRef = useRef(null);
    useEffect(() => {
        plotDataRef.current = plotData;
    });
    const inFlightRef = useRef(false);
    const abortRef = useRef(null);
    // Last LLM call bookkeeping for the throttle gate.
    const lastLlmRef = useRef({ status: null, stats: null, at: 0 });

    useEffect(() => () => abortRef.current?.abort(), []);

    const { baseUrl, model, apiKey, watchdogIntervalMin = 5 } = settings || {};

    useEffect(() => {
        if (!enabled || nSteps < 2) return;
        const status = classifyRun(plotDataRef.current);
        const stats = watchdogRunStats(plotDataRef.current);

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
            apiKey,
            messages: buildWatchdogMessages(stats, status, lastLlmRef.current.status, column),
            temperature: temperatureFor({ baseUrl, model }, 0),
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
    }, [enabled, nSteps, lastValue, baseUrl, model, apiKey, watchdogIntervalMin, column]);

    if (!enabled || nSteps < 2) {
        return { status: 'off', source: 'heuristic', note: null, lastCheckedAt: null, column, allTerms };
    }
    return { ...state, column, allTerms };
};
