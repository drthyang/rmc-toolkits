// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

import { useCallback, useEffect, useRef, useState } from 'react';
import { runCopilotTurn } from './copilot/agentLoop';
import { applyTurnEvent, closeTurn, newAssistantTurn } from './copilot/turnState';
import { buildCopilotMessages, COPILOT_PROMPT_NO_TOOLS } from './prompts/copilot';
import { streamChat } from './provider/client';

// Conversation state for the Copilot page: the finished turns, the turn in
// progress, Stop, and which models proved unable to call tools (so later
// questions skip straight to plain answers). Each question runs one agent
// turn (copilot/agentLoop.js) against the run context and the tools captured
// when it was asked — a Live Data refresh mid-turn changes neither. Opening a
// different run (`runKey`) stops any turn and starts a new conversation, so
// answers about one run never travel as memory into another.
//
// `focus` (the dock's card-focus controller, focus/focus.js) is handed to the
// tools' run() for the view tools; a new question clears the last answer's
// highlight, and the cards of a running analysis breathe until it ends.
// `onTurnEnd(turn)` reports each finished answer (the dock's unread mark).

const now = () => (typeof performance !== 'undefined' ? performance.now() : Date.now());

/**
 * Earlier turns as model memory: each question with its final answer. A
 * question whose turn produced no answer (stopped, failed) is dropped with it,
 * so the memory never holds two user messages in a row or an empty reply.
 */
export const conversationMemory = (turns) => {
    const memory = [];
    turns.forEach((turn, index) => {
        if (turn.role !== 'user') return;
        const reply = turns[index + 1];
        const answer = reply?.role === 'assistant' ? (reply.text || '').trim() : '';
        if (!answer) return;
        memory.push({ role: 'user', content: turn.content }, { role: 'assistant', content: answer });
    });
    return memory;
};

export const useCopilot = ({ settings, context, tools, toolsEnabled = true, runKey = null, focus = null, onTurnEnd = null }) => {
    const [turns, setTurns] = useState([]);
    const [live, setLive] = useState(null);
    const [error, setError] = useState(null);
    const [unsupportedModels, setUnsupportedModels] = useState(() => new Set());
    const abortRef = useRef(null);
    const onTurnEndRef = useRef(onTurnEnd);
    useEffect(() => { onTurnEndRef.current = onTurnEnd; }, [onTurnEnd]);
    // Bumped by clear(): a turn still finishing after it must not reappear.
    const generationRef = useRef(0);

    useEffect(() => () => abortRef.current?.abort(), []);

    const stop = useCallback(() => abortRef.current?.abort(), []);

    const toolsOn = toolsEnabled && tools.length > 0 && !unsupportedModels.has(settings.model);

    const ask = useCallback(async (question) => {
        const text = question.trim();
        if (!text || abortRef.current) return;
        const controller = new AbortController();
        abortRef.current = controller;
        const generation = generationRef.current;
        setError(null);

        // Earlier turns travel as their question and final answer only.
        const history = conversationMemory(turns);
        const viewTools = toolsOn && tools.some((tool) => tool?.view);
        const messages = buildCopilotMessages({ context, history, userText: text, toolsEnabled: toolsOn, viewTools });
        // The agent moves on: the last answer's highlight ends.
        focus?.clear('question');
        const model = settings.model;
        const callModel = (request) => streamChat({
            baseUrl: settings.baseUrl,
            model,
            temperature: settings.temperature,
            apiKey: settings.apiKey,
            ...request
        });

        setTurns((current) => [...current, { role: 'user', content: text }]);
        let turn = newAssistantTurn(now());
        setLive(turn);
        try {
            for await (const event of runCopilotTurn({
                callModel,
                messages,
                tools,
                signal: controller.signal,
                toolsEnabled: toolsOn,
                noToolsSystemPrompt: COPILOT_PROMPT_NO_TOOLS,
                toolContext: focus ? { focus } : null
            })) {
                turn = applyTurnEvent(turn, event, now());
                if (event.type === 'tool_start' && event.cards?.length) focus?.running(event.id, event.cards);
                if (event.type === 'tool_end') focus?.done(event.id);
                if (event.type === 'tools_unsupported') {
                    setUnsupportedModels((current) => new Set(current).add(model));
                }
                setLive(turn);
            }
        } catch (failure) {
            if (failure?.name === 'AbortError') {
                turn = closeTurn(turn, 'stopped');
            } else {
                turn = closeTurn(turn, 'error');
                setError(failure?.message || 'The model request failed');
            }
        } finally {
            abortRef.current = null;
            focus?.doneAll();
            // An empty reply is kept too, so its "The model returned an empty
            // reply" note shows under the question instead of nothing.
            const hasContent = turn.items.some((item) => item.kind === 'tool' || item.text || item.reasoning)
                || turn.stoppedBy === 'empty';
            if (hasContent && generation === generationRef.current) {
                setTurns((current) => [...current, turn]);
                onTurnEndRef.current?.(turn);
            }
            setLive(null);
        }
    }, [turns, context, tools, toolsOn, focus, settings.baseUrl, settings.model, settings.temperature, settings.apiKey]);

    const clear = useCallback(() => {
        generationRef.current += 1;
        abortRef.current?.abort();
        setTurns([]);
        setError(null);
        focus?.clear('new-conversation');
    }, [focus]);

    // Another run opened: stop and forget, as a new conversation.
    const lastRunKeyRef = useRef(runKey);
    useEffect(() => {
        if (lastRunKeyRef.current === runKey) return;
        lastRunKeyRef.current = runKey;
        clear();
    }, [runKey, clear]);

    return {
        turns,
        live,
        busy: Boolean(live),
        error,
        ask,
        stop,
        clear,
        toolsActive: toolsOn,
        toolsUnsupported: toolsEnabled && tools.length > 0 && unsupportedModels.has(settings.model)
    };
};
