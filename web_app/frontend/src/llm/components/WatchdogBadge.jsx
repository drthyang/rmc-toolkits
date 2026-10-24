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
    // The status judges every chi^2 term of the .log and their total when the
    // log names its columns (watchdog/heuristics.js); otherwise ONE fit
    // term's chi^2 (the last .log column), and the badge names that term
    // instead of implying the total.
    const term = !watch.allTerms && watch.column ? `${watch.column}: ` : '';
    let scope = watch.column ? ` — χ² of ${watch.column} only, not the total` : '';
    if (watch.allTerms) scope = ' — every χ² term of the .log and their total: converged only when none still falls beyond its noise over the last tenth';
    const title = watch.note
        ? `${watch.note} — ${source}${scope}`
        : `Convergence watchdog (${source})${scope}`;

    return (
        <span className={`llm-watchdog-badge is-${watch.status}`} title={title} role="status">
            {term}{STATUS_LABELS[watch.status] || watch.status}
            {watch.source === 'llm' && <span className="llm-watchdog-source">AI</span>}
        </span>
    );
};

export default WatchdogBadge;
