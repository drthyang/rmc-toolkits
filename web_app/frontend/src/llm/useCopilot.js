// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

import { useCallback, useEffect, useMemo, useRef, useState } from 'react';
import { describeTurn, runTurn } from './copilot/turn';
import { stoppedByFailure } from './copilot/turnPolicy';
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
// A new configuration of the same run (`configKey`: a Live Data save brought
// a new .rmc6f) keeps the conversation on screen, but its answers were about
// the previous configuration: each question remembers the configuration it
// was asked under, the dock draws a divider where that changes
// (configurationBreaks), and a turn's history holds only the questions asked
// under the current configuration (turnHistory), so the model never reads an
// answer about the previous one as current. A turn running at the save was
// asked about the previous configuration (its run was captured then): the
// divider follows its answer.
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
const EMPTY_DESCRIPTION = Object.freeze({
    offered: Object.freeze({ tools: [], skills: [], askUser: false, view: false }),
    skillsError: null
});

// The configuration a question was asked under (null: the host names none).
const askedUnder = (turn) => turn.configKey ?? null;

/**
 * Where the conversation crosses a new configuration of its run: `before`,
 * the indices of the questions asked under another configuration than the
 * question before them, and `after`, whether the last question was asked
 * under another configuration than `configKey`, the current one. The dock
 * draws a divider before each of `before`, and at the end when `after`.
 */
export const configurationBreaks = (turns, configKey = null) => {
    const before = [];
    let previous;
    turns.forEach((turn, index) => {
        if (turn.role !== 'user') return;
        const key = askedUnder(turn);
        if (previous !== undefined && key !== previous) before.push(index);
        previous = key;
    });
    return { before, after: previous !== undefined && previous !== configKey };
};

/**
 * The conversation's earlier questions as the turn's history: each question
 * with the text of the assistant turn that followed it ('' when none did: a
 * turn that failed before any output is not kept). The turn drops the
 * unanswered ones and compacts long answers. Only the questions asked under
 * `configKey`, the current configuration, since the last change count: an
 * answer about a previous configuration is never sent as memory.
 */
export const turnHistory = (turns, configKey = null) => {
    const { before, after } = configurationBreaks(turns, configKey);
    if (after) return [];
    const current = turns.slice(before.at(-1) ?? 0);
    return current.flatMap((turn, index) => {
        if (turn.role !== 'user') return [];
        const reply = current[index + 1];
        return [{ question: turn.content, answer: reply?.role === 'assistant' ? reply.text || '' : '' }];
    });
};

export const useCopilot = ({
    settings, run, tools, toolsEnabled = true, runKey = null, configKey = null, focus = null, onTurnEnd = null,
    onQuestion = null, skills = null
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

    // What the Copilot can run on this run, with its analyses on. A host that
    // hands the turn what it refuses (a tool named ask_user, say) gets the
    // reason in the error line instead of a broken panel; asking then fails
    // with the same reason.
    const description = useMemo(() => {
        try {
            return describeTurn({ run, tools, askUser, skills, config: { toolsEnabled: true } });
        } catch (failure) {
            return { ...EMPTY_DESCRIPTION, setupError: failure?.message || String(failure) };
        }
    }, [run, tools, askUser, skills]);
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

        setTurns((current) => [...current, { role: 'user', content: text, configKey }]);
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
                // Earlier turns travel as their question and final answer
                // only, and only those about the current configuration.
                history: turnHistory(turns, configKey),
                signal: controller.signal,
                clock: CLOCK,
                toolContext: focus ? { focus } : null,
                config: { toolsEnabled: toolsEnabled && !unsupported }
            })) {
                if (event.type === 'stopped' || event.type === 'failed') {
                    // The turn's ending (turnPolicy.js STOPPED_BY), as the
                    // harness records it.
                    turn = closeTurn(turn, event.stoppedBy);
                    if (event.type === 'failed') setError(event.hint ? `${event.message} ${event.hint}` : event.message);
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
            const stoppedBy = stoppedByFailure(failure);
            turn = closeTurn(turn, stoppedBy);
            if (stoppedBy === 'error') setError(failure?.message || String(failure));
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
    }, [turns, configKey, run, tools, toolsEnabled, unsupported, focus, askUser, skills, settings]);

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

    // Where the conversation crosses a new configuration (the dock's dividers).
    const breaks = useMemo(() => configurationBreaks(turns, configKey), [turns, configKey]);

    return {
        turns,
        breaks,
        live,
        busy: Boolean(live),
        error: error ?? description.setupError ?? null,
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
