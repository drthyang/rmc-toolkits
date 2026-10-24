// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

import { SYSTEM_PROMPT } from './system';

// Message builders for the convergence watchdog. (The Copilot's messages are
// built in prompts/copilot.js.)

export const WATCHDOG_STATUSES = ['improving', 'converged', 'stalled', 'diverging'];

// The watchdog gets a minimal context (recent-window stats, not the full run)
// and must answer in a rigid one-line format so the reply is machine-parseable.
// Stats of the total (`of: 'total'`, heuristics.js watchdogRunStats) are the
// sum of every chi^2 term, with the terms still moving named; otherwise
// `column` names the .log column the statistics come from (e.g. X_ray_(R)1):
// ONE fit term's chi^2, not the total, so the reply must not read it as the
// whole run's convergence.
export const buildWatchdogMessages = (stats, heuristicStatus, prevStatus, column = null) => [
    { role: 'system', content: SYSTEM_PROMPT },
    {
        role: 'user',
        content: 'You are watching a live RMC modeling run. Recent convergence statistics '
            + (stats?.of === 'total'
                ? '(values are ln of the total chi^2, the sum of every chi^2 term of the .log; still_falling and rising name '
                    + 'what still moves beyond its noise over the last tenth; '
                : `(values are ln of chi^2${column ? ` of the fit term ${column} — one .log column, not the total` : ''}, `)
            + `lower is better): ${JSON.stringify(stats)}. `
            + `A simple slope heuristic classifies this as "${heuristicStatus}"`
            + `${prevStatus ? `; the previous assessment was "${prevStatus}"` : ''}. `
            + 'Reply with EXACTLY one line in the form '
            + '`STATUS: improving|converged|stalled|diverging — <one short sentence citing a number>` '
            + 'and nothing else.'
    }
];

// Parse the rigid watchdog reply; returns { status, note } or null when the
// model did not follow the format (callers fall back to the heuristic label).
export const parseWatchdogReply = (text) => {
    const match = /STATUS:\s*(improving|converged|stalled|diverging)\s*[—–-]?\s*(.*)/i.exec(text || '');
    if (!match) return null;
    return { status: match[1].toLowerCase(), note: match[2].trim() || null };
};
