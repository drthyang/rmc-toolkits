// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

import { SYSTEM_PROMPT } from './system';
import { contextToJson } from '../context/runContext';

// Message builders for each assistant feature. Each returns the `messages`
// array for a chat-completions call; the run context always travels inside a
// fenced JSON block so the model can tell data from instructions.

const contextBlock = (context) => `Context for this RMC modeling run:\n\`\`\`json\n${contextToJson(context)}\n\`\`\``;

// Keep only the newest turns so long conversations stay inside small local
// models' context windows; the run context is re-sent every call anyway.
export const CHAT_HISTORY_TURNS = 8;

export const buildChatMessages = (context, history, userText) => [
    { role: 'system', content: SYSTEM_PROMPT },
    { role: 'user', content: contextBlock(context) },
    {
        role: 'assistant',
        content: 'Understood. I will answer questions about this RMC modeling run using only the context above.'
    },
    ...(history || []).slice(-CHAT_HISTORY_TURNS).map(({ role, content }) => ({ role, content })),
    { role: 'user', content: userText }
];

export const WATCHDOG_STATUSES = ['improving', 'converged', 'stalled', 'diverging'];

// The watchdog gets a minimal context (recent-window stats, not the full run)
// and must answer in a rigid one-line format so the reply is machine-parseable.
export const buildWatchdogMessages = (stats, heuristicStatus, prevStatus) => [
    { role: 'system', content: SYSTEM_PROMPT },
    {
        role: 'user',
        content: 'You are watching a live RMC modeling run. Recent convergence statistics '
            + `(values are ln of chi^2, lower is better): ${JSON.stringify(stats)}. `
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
