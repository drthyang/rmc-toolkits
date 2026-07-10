// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

import React from 'react';
import { useLlmSettings } from '../settings';
import { useWatchdog } from '../watchdog/useWatchdog';

// Convergence status chip for the dashboard toolbar. Driven by the heuristics
// in watchdog/heuristics.js (always available) with an LLM-written note when a
// local model is connected. Renders nothing until the user enables the
// watchdog in the assistant's settings.

const STATUS_LABELS = {
    improving: 'Improving',
    converged: 'Converged',
    stalled: 'Stalled',
    diverging: 'Diverging',
    unknown: 'Watching'
};

const WatchdogBadge = ({ rValueFile }) => {
    const settings = useLlmSettings();
    const watch = useWatchdog({ rValueFile, settings });
    if (watch.status === 'off') return null;

    const source = watch.source === 'llm' ? settings.model || 'LLM' : 'heuristic';
    const title = watch.note
        ? `${watch.note} — ${source}`
        : `Convergence watchdog (${source})`;

    return (
        <span className={`llm-watchdog-badge is-${watch.status}`} title={title} role="status">
            {STATUS_LABELS[watch.status] || watch.status}
            {watch.source === 'llm' && <span className="llm-watchdog-source">AI</span>}
        </span>
    );
};

export default WatchdogBadge;
