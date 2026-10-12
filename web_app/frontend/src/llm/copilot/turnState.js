// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// The transcript of one assistant turn, folded from the agent loop's events
// (agentLoop.js). A turn is an ordered list of items — text segments (one per
// model request, with that request's reasoning) and tool steps — plus the
// final verdict. Pure: `applyTurnEvent(turn, event, now)` returns a new turn.

export const newAssistantTurn = (now = 0) => ({
    role: 'assistant',
    items: [],
    text: '',
    outcome: null,
    stoppedBy: null,
    toolsUnsupported: false,
    streaming: true,
    startedAt: now
});

const updateLastText = (items, patch) => {
    const index = items.map((item) => item.kind).lastIndexOf('text');
    if (index < 0) return items;
    const next = [...items];
    next[index] = { ...next[index], ...patch(next[index]) };
    return next;
};

const closeTexts = (items) => items.map((item) => (item.kind === 'text' && item.streaming ? { ...item, streaming: false } : item));

export const applyTurnEvent = (turn, event, now = 0) => {
    switch (event.type) {
        case 'round':
            return {
                ...turn,
                items: [...closeTexts(turn.items), { kind: 'text', text: '', reasoning: '', reasoningMs: null, streaming: true, startedAt: now }]
            };
        case 'text':
            return {
                ...turn,
                items: updateLastText(turn.items, (item) => ({
                    text: item.text + event.delta,
                    // The first answer token freezes how long the model thought.
                    reasoningMs: !item.text && item.reasoning && item.reasoningMs == null ? now - item.startedAt : item.reasoningMs
                }))
            };
        case 'retract':
            // A call written as text: the text was the call, not an answer.
            return { ...turn, items: updateLastText(turn.items, () => ({ text: '' })) };
        case 'reasoning':
            return { ...turn, items: updateLastText(turn.items, (item) => ({ reasoning: item.reasoning + event.delta })) };
        case 'tool_start':
            return {
                ...turn,
                items: [...closeTexts(turn.items), {
                    kind: 'tool',
                    id: event.id,
                    name: event.name,
                    label: event.label,
                    args: event.args,
                    page: event.page,
                    interactive: Boolean(event.interactive),
                    app: Boolean(event.app),
                    reading: Boolean(event.reading),
                    status: 'running',
                    data: null,
                    checks: [],
                    error: null,
                    ms: null,
                    cached: false,
                    pageParams: null
                }]
            };
        case 'tool_end':
            return {
                ...turn,
                items: turn.items.map((item) => (item.kind === 'tool' && item.id === event.id ? {
                    ...item,
                    status: event.status,
                    args: event.args ?? item.args,
                    data: event.data,
                    checks: event.checks || [],
                    error: event.error,
                    ms: event.ms,
                    cached: Boolean(event.cached),
                    pageParams: event.pageParams ?? null
                } : item))
            };
        case 'tools_unsupported':
            return { ...turn, toolsUnsupported: true };
        case 'done': {
            // The verdict line is shown as the chip, not repeated in the text.
            const items = event.outcome
                ? updateLastText(turn.items, () => ({ text: event.outcome.body }))
                : turn.items;
            return {
                ...turn,
                items: closeTexts(items),
                text: event.text,
                outcome: event.outcome,
                stoppedBy: event.stoppedBy,
                cut: Boolean(event.cut),
                reason: event.reason || null,
                audit: Array.isArray(event.audit) ? event.audit : [],
                streaming: false
            };
        }
        default:
            return turn;
    }
};

/**
 * A turn interrupted by Stop or an error: close it and mark running steps
 * stopped. Its last partial text becomes the turn's memory, marked as
 * unfinished so a later question does not read it as a complete answer.
 */
export const closeTurn = (turn, stoppedBy) => {
    const partial = turn.text || [...turn.items].reverse().find((item) => item.kind === 'text' && item.text)?.text || '';
    return {
        ...turn,
        items: closeTexts(turn.items).map((item) => (item.kind === 'tool' && item.status === 'running'
            ? { ...item, status: 'stopped' }
            : item)),
        text: partial.trim() ? `${partial.trim()}\n\n(Interrupted before the answer was finished.)` : '',
        stoppedBy,
        streaming: false
    };
};

/**
 * "3 analyses · 1 warning · 1 failed" for a turn's footer, or '' without
 * steps. A clarifying question (an interactive step) and reading a skill are
 * not analyses.
 */
export const stepSummary = (items) => {
    const steps = items.filter((item) => item.kind === 'tool' && !item.interactive && !item.reading);
    if (!steps.length) return '';
    const count = (status) => steps.filter((item) => item.status === status).length;
    const parts = [`${steps.length} ${steps.length === 1 ? 'analysis' : 'analyses'}`];
    const warn = count('warn');
    const failed = count('failed');
    if (warn) parts.push(`${warn} ${warn === 1 ? 'warning' : 'warnings'}`);
    if (failed) parts.push(`${failed} failed`);
    return parts.join(' · ');
};
