// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

import { useCallback, useEffect, useMemo, useRef, useState } from 'react';
import { describeTurn, runTurn } from './copilot/turn';
import { applyTurnEvent, closeTurn, newAssistantTurn } from './copilot/turnState';
import { providerModel } from './provider/model';

// The dock's adapter of the turn module (copilot/turn.js): the conversation
// state for the Copilot page — the finished turns, the turn in progress,
// Stop, and which models proved unable to call tools (so later questions
// skip straight to plain answers). Each question runs one turn against the
// run props and the tools captured when it was asked — a Live Data refresh
// mid-turn changes neither. The turn builds everything else (the context, the
// prompt, ask_user and read_skill, the checks); this hook turns its events
// into React state. Opening a different run (`runKey`) stops any turn and
// starts a new conversation, so answers about one run never travel as memory
// into another.
//
// `focus` (the dock's card-focus controller, focus/focus.js) is handed to the
// tools' run() for the view tools; a new question clears the last answer's
// highlight, and the cards of a running analysis breathe until it ends.
// `onTurnEnd(turn)` reports each finished answer (the dock's unread mark).
//
// ask_user (copilot/askUser.js) reaches the user through `askUser`, the
// turn's ask-user port: it holds the question the card shows (`question`)
// until `answerQuestion(reply)` settles it, and rejects it with an AbortError
// when the turn is stopped, the conversation restarts or the panel unmounts.
// `onQuestion(request)` reports a new question (the dock's unread mark).
//
// Every model request goes through providerModel (provider/model.js) and so
// through provider/contextGuard.js, which never sends a request the model
// would cut to fit its context window (it sizes Ollama's context, or refuses
// the turn with a message that says why); `session` keeps the context asked
// per model across the conversation.
//
// `offered` and `skillsError` describe what the Copilot can run on this run
// with its analyses on (describeTurn): the composer's status line and the
// tool list. Whether a question uses them is the analyses switch
// (`toolsEnabled`) and the models that refused tools.

const now = () => (typeof performance !== 'undefined' ? performance.now() : Date.now());
const CLOCK = Object.freeze({ now });
const abortError = () => new DOMException('Aborted', 'AbortError');

/**
 * The conversation's earlier questions as the turn's history: each question
 * with the text of the assistant turn that followed it ('' when none did: a
 * turn that failed before any output is not kept). The turn drops the
 * unanswered ones and compacts long answers.
 */
export const turnHistory = (turns) => turns.flatMap((turn, index) => {
    if (turn.role !== 'user') return [];
    const reply = turns[index + 1];
    return [{ question: turn.content, answer: reply?.role === 'assistant' ? reply.text || '' : '' }];
});

export const useCopilot = ({
    settings, run, tools, toolsEnabled = true, runKey = null, focus = null, onTurnEnd = null, onQuestion = null,
    skills = null
}) => {
    const [turns, setTurns] = useState([]);
    const [live, setLive] = useState(null);
    const [error, setError] = useState(null);
    const [unsupportedModels, setUnsupportedModels] = useState(() => new Set());
    // The question waiting for the user: { request, reply(answer) }.
    const [question, setQuestion] = useState(null);
    const questionRef = useRef(null);
    const abortRef = useRef(null);
    const onTurnEndRef = useRef(onTurnEnd);
    useEffect(() => { onTurnEndRef.current = onTurnEnd; }, [onTurnEnd]);
    const onQuestionRef = useRef(onQuestion);
    useEffect(() => { onQuestionRef.current = onQuestion; }, [onQuestion]);
    // Bumped by clear(): a turn still finishing after it must not reappear.
    const generationRef = useRef(0);
    // The context asked for per model, kept across this conversation's turns
    // so a model is not reloaded for every question (contextGuard.js).
    const contextSessionRef = useRef({});

    useEffect(() => () => abortRef.current?.abort(), []);

    const stop = useCallback(() => abortRef.current?.abort(), []);

    // ask_user's way to the user: show the question, wait for the reply.
    const askUser = useCallback((request, { signal } = {}) => new Promise((resolve, reject) => {
        if (signal?.aborted) {
            reject(abortError());
            return;
        }
        const pending = { request };
        const settle = () => {
            signal?.removeEventListener('abort', onAbort);
            if (questionRef.current === pending) {
                questionRef.current = null;
                setQuestion(null);
            }
        };
        const onAbort = () => {
            settle();
            reject(abortError());
        };
        pending.reply = (reply) => {
            settle();
            resolve(reply);
        };
        signal?.addEventListener('abort', onAbort, { once: true });
        questionRef.current = pending;
        setQuestion(pending);
        onQuestionRef.current?.(request);
    }), []);

    /** The user's reply to the waiting question: { answer, source, index }. */
    const answerQuestion = useCallback((reply) => questionRef.current?.reply(reply), []);

    // What the Copilot can run on this run, with its analyses on.
    const description = useMemo(
        () => describeTurn({ run, tools, askUser, skills, config: { toolsEnabled: true } }),
        [run, tools, askUser, skills]
    );
    const offersTools = description.offered.tools.length > 0;
    const unsupported = unsupportedModels.has(settings.model);
    const toolsOn = toolsEnabled && offersTools && !unsupported;

    const ask = useCallback(async (questionText) => {
        const text = questionText.trim();
        if (!text || abortRef.current) return;
        const controller = new AbortController();
        abortRef.current = controller;
        const generation = generationRef.current;
        setError(null);

        // The agent moves on: the last answer's highlight ends.
        focus?.clear('question');
        const model = settings.model;

        setTurns((current) => [...current, { role: 'user', content: text }]);
        let turn = newAssistantTurn(now());
        setLive(turn);
        try {
            for await (const event of runTurn({
                question: text,
                run,
                model: providerModel(settings, { session: contextSessionRef.current }),
                tools,
                askUser,
                skills,
                // Earlier turns travel as their question and final answer only.
                history: turnHistory(turns),
                signal: controller.signal,
                clock: CLOCK,
                toolContext: focus ? { focus } : null,
                config: { toolsEnabled: toolsEnabled && !unsupported }
            })) {
                if (event.type === 'stopped') {
                    turn = closeTurn(turn, 'stopped');
                } else if (event.type === 'failed') {
                    turn = closeTurn(turn, 'error');
                    setError(event.hint ? `${event.message} ${event.hint}` : event.message);
                } else {
                    turn = applyTurnEvent(turn, event, event.at);
                    if (event.type === 'tool_start' && event.cards?.length) focus?.running(event.id, event.cards);
                    if (event.type === 'tool_end') focus?.done(event.id);
                    if (event.type === 'tools_unsupported') {
                        setUnsupportedModels((current) => new Set(current).add(model));
                    }
                }
                setLive(turn);
            }
        } catch (failure) {
            // The turn never throws once started: this is the dock's own
            // failure, shown as any other.
            turn = closeTurn(turn, failure?.name === 'AbortError' ? 'stopped' : 'error');
            if (failure?.name !== 'AbortError') setError(failure?.message || String(failure));
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
    }, [turns, run, tools, toolsEnabled, unsupported, focus, askUser, skills, settings]);

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
        question: question?.request || null,
        answerQuestion,
        offered: description.offered,
        skillsError: description.skillsError,
        toolsActive: toolsOn,
        toolsUnsupported: toolsEnabled && offersTools && unsupported
    };
};
