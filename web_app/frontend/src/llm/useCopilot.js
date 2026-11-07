// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

import { useCallback, useEffect, useMemo, useRef, useState } from 'react';
import { compactAnswer } from './copilot/robustness';
import { MEMORY_ANSWER_CHARS, describeTurn, historyWindow, requestPreview, runTurn } from './copilot/turn';
import { stoppedByFailure } from './copilot/turnPolicy';
import { applyTurnEvent, closeTurn, newAssistantTurn } from './copilot/turnState';
import { CHARS_PER_TOKEN, estimatePromptTokens } from './provider/contextGuard';
import { contextUsage, meterLimits, probeLimits, requestCost } from './provider/contextMeter';
import { providerModel } from './provider/model';
import { temperatureFor } from './provider/sampling';

// Summarize and continue: each earlier answer is re-sent as at most this
// many characters (compactAnswer: its sentences with numbers, its Outcome
// line). No model call writes the summary, so it cannot add or change a
// number.
export const SUMMARY_ANSWER_CHARS = 600;

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
// A change tool (copilot/effects.js) asks first through `approve`, the
// turn's approve port: with `autoApprove` off (Ask first, the default) it
// holds the change the card shows (`approval`) until `decideChange(approved)`
// settles it ({ approved, auto: false }), and rejects with an AbortError as a
// question does; with `autoApprove` on (Auto, which the dock resets on every
// reload) it approves at once ({ approved: true, auto: true }).
// `onApproval(request)` reports a change waiting (the unread mark).
// `undoChange(step)` undoes an applied change through its tool's `undo`; how
// it went is `undos[change id]` ({ state: 'running' | 'done' | 'failed',
// reason }), keyed by the applied change's id (unique in the app, where a
// step's id may repeat in a later turn), cleared with the conversation.
//
// Every model request goes through providerModel (provider/model.js) and so
// through provider/contextGuard.js, which never sends a request the model
// would cut to fit its context window (it sizes Ollama's context, or refuses
// the turn with a message that says why); `session` keeps the context asked
// per model across the conversation. A model that takes no temperature
// (provider/sampling.js: Claude) is given `temperature: null`, so none of its
// requests carries one, the warmer retry included (copilotModel).
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

/**
 * The dock's model port from the user's settings: providerModel with the
 * settings' temperature, or null (none sent, the warmer retry included) for
 * a model that takes none (provider/sampling.js takesNoTemperature).
 */
export const copilotModel = (settings, options) => providerModel(
    { ...settings, temperature: temperatureFor(settings) },
    options
);

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
export const turnHistory = (turns, configKey = null, { summarizedAt = null } = {}) => {
    const { before, after } = configurationBreaks(turns, configKey);
    if (after) return [];
    const start = before.at(-1) ?? 0;
    const current = turns.slice(start);
    return current.flatMap((turn, index) => {
        if (turn.role !== 'user') return [];
        const reply = current[index + 1];
        const answer = reply?.role === 'assistant' ? reply.text || '' : '';
        // Summarize and continue: the answers before that point go as their
        // sentences with numbers (and their Outcome line), shorter than the
        // turn's own MEMORY_ANSWER_CHARS.
        const summarized = Number.isInteger(summarizedAt) && start + index < summarizedAt;
        return [{ question: turn.content, answer: summarized && answer ? compactAnswer(answer, SUMMARY_ANSWER_CHARS) : answer }];
    });
};

// The output tokens of a reply, estimated as the guard estimates text.
const replyTokens = (chars) => Math.ceil(chars / CHARS_PER_TOKEN);

// A request's tokens as the guard counts them: the plan's own estimate on
// Ollama and LM Studio (the same estimatePromptTokens), else counted here.
const requestEstimate = (sent) => (Number.isFinite(sent?.plan?.estimate)
    ? sent.plan.estimate
    : estimatePromptTokens({ messages: sent?.messages || [], tools: sent?.tools ?? null }));

// The larger of two requests (by tokens).
const larger = (a, b) => (!a ? b : (!b ? a : (b.estimate > a.estimate ? b : a)));

export const useCopilot = ({
    settings, run, tools, toolsEnabled = true, runKey = null, configKey = null, focus = null, onTurnEnd = null,
    onQuestion = null, onApproval = null, skills = null, autoApprove = false, connected = false
}) => {
    const [turns, setTurns] = useState([]);
    const [live, setLive] = useState(null);
    const [error, setError] = useState(null);
    const [unsupportedModels, setUnsupportedModels] = useState(() => new Set());
    // The question waiting for the user: { request, reply(answer) }.
    const [question, setQuestion] = useState(null);
    const questionRef = useRef(null);
    // The change waiting for the user's approval: { request, decide(approved) }.
    const [approval, setApproval] = useState(null);
    const approvalRef = useRef(null);
    // Read when a change asks, so a switch to Auto applies from the next change on.
    const autoApproveRef = useRef(autoApprove);
    autoApproveRef.current = autoApprove;
    // How each Undo went, by the applied change's id.
    const [undos, setUndos] = useState({});
    const toolsRef = useRef(tools);
    toolsRef.current = tools;
    const abortRef = useRef(null);
    const onTurnEndRef = useRef(onTurnEnd);
    useEffect(() => { onTurnEndRef.current = onTurnEnd; }, [onTurnEnd]);
    const onQuestionRef = useRef(onQuestion);
    useEffect(() => { onQuestionRef.current = onQuestion; }, [onQuestion]);
    const onApprovalRef = useRef(onApproval);
    useEffect(() => { onApprovalRef.current = onApproval; }, [onApproval]);
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

    // A change tool's way to the user: Auto approves at once; Ask first shows
    // the card and waits for Apply or Decline.
    const approve = useCallback((request, { signal } = {}) => new Promise((resolve, reject) => {
        if (signal?.aborted) {
            reject(abortError());
            return;
        }
        if (autoApproveRef.current) {
            resolve({ approved: true, auto: true });
            return;
        }
        const pending = { request };
        const settle = () => {
            signal?.removeEventListener('abort', onAbort);
            if (approvalRef.current === pending) {
                approvalRef.current = null;
                setApproval(null);
            }
        };
        const onAbort = () => {
            settle();
            reject(abortError());
        };
        pending.decide = (approved) => {
            settle();
            resolve({ approved: Boolean(approved), auto: false });
        };
        signal?.addEventListener('abort', onAbort, { once: true });
        approvalRef.current = pending;
        setApproval(pending);
        onApprovalRef.current?.(request);
    }), []);

    /** The user's decision on the waiting change: Apply (true) or Decline (false). */
    const decideChange = useCallback((approved) => approvalRef.current?.decide(approved), []);

    /** Undo an applied change (a change step's record) through its tool's `undo`. */
    const undoChange = useCallback(async (step) => {
        const tool = (toolsRef.current || []).find((item) => item?.name === step?.name);
        if (!step?.change?.id || typeof tool?.undo !== 'function') return;
        const key = step.change.id;
        setUndos((current) => ({ ...current, [key]: { state: 'running' } }));
        let outcome;
        try {
            outcome = await tool.undo(step.change);
        } catch (failure) {
            outcome = { ok: false, reason: failure?.message || String(failure) };
        }
        setUndos((current) => ({
            ...current,
            [key]: outcome?.ok ? { state: 'done' } : { state: 'failed', reason: outcome?.reason || 'it could not be undone' }
        }));
    }, []);

    // What the Copilot can run on this run, with its analyses on. A host that
    // hands the turn what it refuses (a tool named ask_user, say) gets the
    // reason in the error line instead of a broken panel; asking then fails
    // with the same reason.
    const description = useMemo(() => {
        try {
            return describeTurn({ run, tools, askUser, approve, skills, config: { toolsEnabled: true } });
        } catch (failure) {
            return { ...EMPTY_DESCRIPTION, setupError: failure?.message || String(failure) };
        }
    }, [run, tools, askUser, approve, skills]);
    const offersTools = description.offered.tools.length > 0;
    const unsupported = unsupportedModels.has(settings.model);
    const toolsOn = toolsEnabled && offersTools && !unsupported;

    // The context meter (provider/contextMeter.js): the guard's numbers for
    // the next question's first request, the largest request of the last
    // answer, the running answer's latest request, or a request the guard
    // refused. `limitsPlan`: the guard's latest plan for this model (a
    // request's, or probeLimits' before the first one).
    const [summarizedAt, setSummarizedAt] = useState(null);
    const [limitsPlan, setLimitsPlan] = useState(null);
    const [liveRequest, setLiveRequest] = useState(null);
    const [peak, setPeak] = useState(null);
    const [refused, setRefused] = useState(null);
    const [cost, setCost] = useState(0);
    const [probeTick, setProbeTick] = useState(0);
    // What the last question's history window held, for the chat's memory note.
    const memoryRef = useRef({ answered: 0, dropped: 0 });

    // The window of a local server's model, before any request (the same
    // probes and caches as the guard; a cloud model needs none). Asked again
    // after each answer: LM Studio's model may have been reloaded.
    useEffect(() => {
        if (!connected || !settings.baseUrl || !settings.model) return undefined;
        const controller = new AbortController();
        const model = settings.model;
        probeLimits({ baseUrl: settings.baseUrl, model }, { signal: controller.signal })
            .then((plan) => setLimitsPlan({ plan, model }))
            // An unreachable server: the window stays unknown (the meter says so).
            .catch(() => {});
        return () => controller.abort();
    }, [connected, settings.baseUrl, settings.model, probeTick]);

    const history = useMemo(() => turnHistory(turns, configKey, { summarizedAt }), [turns, configKey, summarizedAt]);
    // The next question's first request, as the guard would estimate it.
    const previewTokens = useCallback((items) => {
        try {
            const { messages, tools: sent } = requestPreview({
                run, tools, askUser, approve, skills, history: items, config: { toolsEnabled: toolsOn }
            });
            return estimatePromptTokens({ messages, tools: sent });
        } catch {
            // The setup error is shown already (description.setupError).
            return null;
        }
    }, [run, tools, askUser, approve, skills, toolsOn]);
    const nextTokens = useMemo(() => previewTokens(history), [previewTokens, history]);
    // What Summarize and continue would leave: offered only when it helps.
    const summarizedTokens = useMemo(
        () => (turns.length ? previewTokens(turnHistory(turns, configKey, { summarizedAt: turns.length })) : null),
        [previewTokens, turns, configKey]
    );

    const ask = useCallback(async (questionText) => {
        const text = questionText.trim();
        if (!text || abortRef.current) return;
        const controller = new AbortController();
        abortRef.current = controller;
        const generation = generationRef.current;
        setError(null);
        setRefused(null);

        // The agent moves on: the last answer's highlight ends.
        focus?.clear('question');
        const model = settings.model;

        // Never truncated silently: what this question's memory leaves out
        // that the last one sent (the oldest questions dropped, the last
        // answer shortened) goes into the chat as the answer's memory note.
        const window = historyWindow(history);
        const known = window.answered >= memoryRef.current.answered ? memoryRef.current : { answered: 0, dropped: 0 };
        const lastAnswer = history.filter((entry) => entry.answer.trim() && entry.question.trim()).at(-1)?.answer.trim() || '';
        const memory = {
            dropped: Math.max(0, window.dropped - known.dropped),
            kept: window.sent,
            shortened: window.answered > known.answered && lastAnswer.length > MEMORY_ANSWER_CHARS
        };
        memoryRef.current = window;

        setTurns((current) => [...current, { role: 'user', content: text, configKey }]);
        let turn = { ...newAssistantTurn(now()), memory: memory.dropped || memory.shortened ? memory : null };
        setLive(turn);
        // This answer's requests: the largest, and its estimated cost.
        let turnPeak = null;
        let spent = 0;
        let replyChars = 0;
        try {
            for await (const event of runTurn({
                question: text,
                run,
                model: copilotModel(settings, { session: contextSessionRef.current }),
                tools,
                askUser,
                approve,
                skills,
                // Earlier turns travel as their question and final answer
                // only, and only those about the current configuration.
                history,
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
                    if (event.error?.name === 'ContextLimitError') {
                        setRefused({ needed: event.error.needed, available: event.error.available });
                    }
                } else {
                    turn = applyTurnEvent(turn, event, event.at);
                    if (event.type === 'tool_start' && event.cards?.length) focus?.running(event.id, event.cards);
                    if (event.type === 'tool_end') focus?.done(event.id);
                    if (event.type === 'tools_unsupported') {
                        setUnsupportedModels((current) => new Set(current).add(model));
                    }
                    if (event.type === 'request') {
                        const request = { estimate: requestEstimate(event.sent) };
                        turnPeak = larger(turnPeak, request);
                        setLiveRequest(request);
                        if (event.sent?.plan) setLimitsPlan({ plan: event.sent.plan, model });
                        spent += requestCost({ model, inputTokens: request.estimate }) ?? 0;
                    }
                    if (event.type === 'text' || event.type === 'reasoning') replyChars += event.delta?.length || 0;
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
            if (generation === generationRef.current) {
                setLiveRequest(null);
                setPeak(turnPeak);
                spent += requestCost({ model, inputTokens: 0, outputTokens: replyTokens(replyChars) }) ?? 0;
                if (spent) setCost((current) => current + spent);
                setProbeTick((tick) => tick + 1);
            }
        }
    }, [configKey, history, run, tools, toolsEnabled, unsupported, focus, askUser, approve, skills, settings]);

    const clear = useCallback(() => {
        generationRef.current += 1;
        abortRef.current?.abort();
        setTurns([]);
        setUndos({});
        setError(null);
        setSummarizedAt(null);
        setLiveRequest(null);
        setPeak(null);
        setRefused(null);
        setCost(0);
        memoryRef.current = { answered: 0, dropped: 0 };
        focus?.clear('new-conversation');
    }, [focus]);

    // Summarize and continue: the answers so far go to the model as their
    // sentences with numbers from the next question on (turnHistory), and
    // the chat shows where (summarizedAt).
    const summarize = useCallback(() => {
        if (abortRef.current || !turns.length) return;
        setSummarizedAt(turns.length);
        setPeak(null);
        setRefused(null);
        setError(null);
    }, [turns.length]);

    // Another run opened: stop and forget, as a new conversation.
    const lastRunKeyRef = useRef(runKey);
    useEffect(() => {
        if (lastRunKeyRef.current === runKey) return;
        lastRunKeyRef.current = runKey;
        clear();
    }, [runKey, clear]);

    // Where the conversation crosses a new configuration (the dock's dividers).
    const breaks = useMemo(() => configurationBreaks(turns, configKey), [turns, configKey]);

    const meter = useMemo(() => {
        const plan = limitsPlan?.model === settings.model ? limitsPlan.plan : null;
        const limits = meterLimits(plan, { model: settings.model });
        const usage = (request) => (request && Number.isFinite(request.estimate) ? contextUsage({ estimate: request.estimate, limits }) : null);
        const next = usage(nextTokens == null ? null : { estimate: nextTokens });
        const last = usage(peak);
        const running = usage(liveRequest);
        let shown = running || (last && next ? (last.used > next.used ? last : next) : (next || last));
        if (refused && Number.isFinite(refused.needed)) {
            // The guard refused: the request it would have needed, against
            // the window it had (100 % and more).
            const window = Number.isFinite(refused.available) && refused.available > 0 ? refused.available : limits.window;
            const share = window ? refused.needed / window : null;
            shown = { prompt: refused.needed - limits.reserve, reserve: limits.reserve, used: refused.needed, window, share, level: 'red', refused: true };
        }
        if (!shown) return null;
        return {
            ...shown,
            next,
            last,
            running,
            limits,
            model: settings.model,
            cost: requestCost({ model: settings.model, inputTokens: 0 }) == null ? null : cost,
            canSummarize: !running && nextTokens != null && summarizedTokens != null && summarizedTokens < nextTokens,
            summarized: summarizedAt != null
        };
    }, [limitsPlan, settings.model, nextTokens, peak, liveRequest, refused, cost, summarizedTokens, summarizedAt]);

    return {
        turns,
        breaks,
        live,
        busy: Boolean(live),
        error: error ?? description.setupError ?? null,
        ask,
        stop,
        clear,
        summarize,
        summarizedAt,
        meter,
        question: question?.request || null,
        answerQuestion,
        approval: approval?.request || null,
        decideChange,
        undos,
        undoChange,
        offered: description.offered,
        skillsError: description.skillsError,
        toolsActive: toolsOn,
        toolsUnsupported: toolsEnabled && offersTools && unsupported
    };
};
