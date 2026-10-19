// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

import { toolBudgetNote } from '../prompts/copilot';

// The turn policy: what ends a turn (architecture candidate #10). One module
// owns every counter that decides when a turn ends — the rounds of tools, the
// turn's time, the rounds that added nothing, the one reminder and the warmer
// retry, the last-request note, the answer the evidence gate withdrew — and the
// order of the steps, behind one event-driven interface:
//
//   createTurnPolicy({ limits, clock, tools, judge }) → { next(event) → step }
//   STOPPED_BY               the one set of endings: whether each cut the analyses
//                            short, and the note the dock shows under the turn
//   stoppedByFailure(error)  the ending of a turn that threw: 'stopped' or 'error'
//
// The agent loop (agentLoop.js) is its one caller and its seam: the loop
// streams the model, runs the calls and keeps the messages; it tells the
// policy what happened (an event) and does what the step says. So the loop
// cannot ask in the wrong order, add the last-request note twice or end a turn
// twice, and the endings are decided in one place: the dock (CopilotView's
// note, useCopilot), the computed Outcome (verification/verdict.js: cutShort)
// and the validation harness (record.stopped_by) all read STOPPED_BY's keys.
//
// Configuration:
//   limits   { maxRounds, maxTurnMs, maxNoProgressRounds }: numbers, Infinity
//            allowed (the loop's DEFAULT_LIMITS fill them); other keys are not
//            read. The send-backs per turn are the verification's own cap.
//   clock    { now() } in ms. The turn starts when its policy is created: a
//            question the Workbench asks before the first request is a wait,
//            reported by `start`.
//   tools    whether the first request carries tools.
//   judge    the verification's judge (verification/): (text, { canSendBack,
//            rounds, stoppedBy, cut }) → { sendBack } | { done } | null. The
//            policy decides WHEN — whether a send-back is allowed, the stoppedBy
//            an answer is judged under, which answer stands — and the
//            verification WHAT.
//
// Events (the loop's), in order: start, then for each request a reply, which
// is followed by a round when it had calls, or no_tools when the model refused
// its tools:
//   { kind: 'start', waitedMs }                   the turn begins; waitedMs: the
//                                                 Workbench's question before it
//   { kind: 'reply', text, calls, looped, cut }   a request's reply: its text
//                                                 (none when looped), the calls
//                                                 it holds (an array), whether it
//                                                 was stopped for repeating itself,
//                                                 and whether it reached the
//                                                 length limit
//   { kind: 'round', progressed, waitedMs }       a round of calls ran: whether it
//                                                 brought anything new, and the
//                                                 time it waited for the user
//   { kind: 'no_tools' }                          the first request's model cannot
//                                                 take tools (the port's
//                                                 ToolsUnsupportedError)
// Steps:
//   { next: 'ask', round, why, toolChoice, retry, notes, sendBack }
//       make a model request. why: 'first' | 'tools' | 'reminder' | 'gate' |
//       'final' | 'no_tools'; toolChoice: 'auto', 'none' (the request that asks
//       for the answer) or undefined (no tools); retry: 'warmer' | null (the
//       model port's); notes: the messages to add (addNote), in order; sendBack:
//       the verification's send-back for the reply just judged ({ id, what,
//       missing, notice, nudge }: the loop withdraws the reply and sends the
//       nudge) or null
//   { next: 'run' }                               run the reply's calls, then
//                                                 report the round
//   { next: 'end', done, say }                    yield `say` as text when it is
//                                                 not null, then `done`; the turn
//                                                 is over
// `done` has one shape on every path, built only here:
//   { type: 'done', text, rawText, outcome, rounds, stoppedBy, cut, reason, audit, references, gates }
//   - an answer (the verification's done): reason null; rawText the answer as
//     the model wrote it (the text itself when no reference changed it);
//   - the request that asked for the answer came back empty, with nothing
//     standing: a fixed answer (`say`), Outcome 'not achieved' with its note and
//     body only, rawText = text, references null, gates [];
//   - an empty reply, after its reminder, before the last request: text '',
//     rawText '', outcome null, stoppedBy 'empty', reason (why it was empty),
//     references null, gates [].
//
// Invariants:
//   - toolChoice is undefined exactly when tools are off;
//   - an ask step's notes hold 0–2 messages, in order (a reminder, then the
//     last-request note); the last-request note is added once; the request
//     that becomes the last gets why 'final', unless it is the first;
//   - one reminder a turn; after a looping reply, every later request is
//     asked warmer;
//   - turn time = now − start − the waits reported, read at most once per
//     event; round 0 is never out of time;
//   - in a reply without calls the order is: the reminder, then the answer the
//     gate withdrew stands, then the empty endings, then the judge;
//   - judge gets only non-empty text; the withdrawn answer with canSendBack false;
//   - canSendBack = tools on ∧ this request is not the last ∧ turn time left
//     (and the verification's own cap);
//   - exactly one 'end' a turn.
// Errors: a TypeError for bad configuration, or for an event out of order or
// malformed (calls on a request that allowed none, no_tools after tools ran),
// a send-back the policy did not allow, or a judge that answers otherwise than
// { sendBack }, { done } or null; judge's own errors propagate (the turn ends
// as `failed`); judge returning null gives the empty ending. Pure, no React.

/** The endings a turn can have: every `stoppedBy` the dock and the harness read. */
export const STOPPED_BY = Object.freeze({
    answer: Object.freeze({ cutShort: false, note: null }),
    // Also a turn out of time: no new round of tools after maxTurnMs.
    max_rounds: Object.freeze({ cutShort: true, note: 'Answered after using the analysis budget for one question.' }),
    no_progress: Object.freeze({ cutShort: true, note: 'Answered when repeated analyses stopped adding anything new.' }),
    empty: Object.freeze({ cutShort: false, note: 'The model returned an empty reply.' }),
    stopped: Object.freeze({ cutShort: false, note: 'Stopped.' }),
    // A stall (ModelTimeoutError) too.
    error: Object.freeze({ cutShort: false, note: 'Stopped by an error.' })
});

/** The ending of a turn that threw: Stop (an AbortError) or any other error. */
export const stoppedByFailure = (error) => (error?.name === 'AbortError' ? 'stopped' : 'error');

const EMPTY_FINAL = 'The analysis budget for this question was used before the model wrote an answer. '
    + 'The analyses it ran are listed above.';
const EMPTY_STALLED = 'The model kept repeating analyses that added nothing new, and wrote no answer. '
    + 'The analyses it ran are listed above.';

// When a reply holds neither an answer nor an analysis call.
const emptyNudge = (why) => ({
    role: 'user',
    content: `Your last reply ${why === 'loop'
        ? 'kept repeating the same passage, so it was stopped'
        : `had no answer and no analysis call${why === 'cut' ? ' (it reached the length limit while reasoning)' : ''}`}. `
        + 'Answer the question now from the context and the results above, briefly, or call one analysis if '
        + 'one is still missing. End with the Outcome line.'
});

const EMPTY_REASONS = {
    loop: 'the model kept repeating the same passage, even after a reminder',
    cut: 'the model spent its reply on reasoning and reached the length limit',
    empty: 'the model wrote no answer, even after a reminder'
};

// When the last rounds added nothing new.
const noProgressNote = (rounds) => ({
    role: 'user',
    content: `Your last ${rounds} steps repeated earlier analyses or failed, so they added nothing new. `
        + 'Do not call them again. Answer now from the results above, say what is still missing, or ask the user '
        + 'with ask_user if a choice blocks you. End with the Outcome line.'
});

const LIMITS = ['maxRounds', 'maxTurnMs', 'maxNoProgressRounds'];

// The events each phase takes: one start, then a request's reply (or the
// model's refusal of its tools), and after a reply with calls their round.
const EXPECTED = Object.freeze({
    start: ['start'],
    asked: ['reply', 'no_tools'],
    running: ['round'],
    ended: []
});

const isObject = (value) => value !== null && typeof value === 'object' && !Array.isArray(value);
const invalid = (message) => {
    throw new TypeError(`turn policy: ${message}`);
};
const isLimit = (value) => typeof value === 'number' && !Number.isNaN(value);
const checkWait = (waitedMs) => {
    if (!Number.isFinite(waitedMs)) invalid(`waitedMs must be a finite number of ms, got ${String(waitedMs)}`);
};

/** One turn's policy (see the header): `next(event)` → the next step. */
export const createTurnPolicy = ({ limits, clock, tools, judge } = {}) => {
    if (!isObject(limits)) invalid('limits must be { maxRounds, maxTurnMs, maxNoProgressRounds }');
    LIMITS.forEach((key) => {
        if (!isLimit(limits[key])) invalid(`limits.${key} must be a number (Infinity allowed), got ${String(limits[key])}`);
    });
    if (!clock || typeof clock.now !== 'function') invalid('clock must have now()');
    if (typeof tools !== 'boolean') invalid('tools must be a boolean: whether the first request carries tools');
    if (typeof judge !== 'function') invalid('judge must be a function: the verification\'s judge');
    const { maxRounds, maxTurnMs, maxNoProgressRounds } = limits;

    const startedAt = clock.now();
    let phase = 'start';
    let toolsOn = tools;
    // Rounds of tool calls run.
    let round = 0;
    // The user's answers (Workbench questions, ask_user), off the turn's time.
    let waited = 0;
    let noProgressRounds = 0;
    let stalled = false;
    let reminded = false;
    let warmer = false;
    let finalNoted = false;
    // Whether the request in flight is the last: it asks for the answer.
    let final = false;
    // The answer the evidence gate sent back last: it stands if the model
    // then writes nothing.
    let withdrawn = null;

    // The next request, from the event's one reading of the turn's time.
    const ask = (reason, outOfTime, { notes = [], sendBack = null } = {}) => {
        final = toolsOn && (round >= maxRounds || stalled || outOfTime());
        let why = reason;
        const carried = [...notes];
        // The last-request note is added once (a reminder can repeat that request).
        if (final && !finalNoted) {
            finalNoted = true;
            carried.push(stalled ? noProgressNote(noProgressRounds) : toolBudgetNote(round));
            // The request that asks for the answer (the first request stays 'first').
            if (why !== 'first') why = 'final';
        }
        phase = 'asked';
        return {
            next: 'ask',
            round,
            why,
            toolChoice: toolsOn ? (final ? 'none' : 'auto') : undefined,
            retry: warmer ? 'warmer' : null,
            notes: carried,
            sendBack
        };
    };

    const end = (done, say = null) => {
        phase = 'ended';
        return { next: 'end', done, say };
    };

    const doneEvent = ({ text, rawText = text, outcome = null, stoppedBy, cut, reason = null, audit = [], references = null, gates = [] }) => ({
        type: 'done', text, rawText, outcome, rounds: round, stoppedBy, cut, reason, audit, references, gates
    });

    // No answer and none standing: on the last request, a fixed answer that
    // says so; before it, the empty ending with why the reply was empty.
    const emptyEnding = (why, stoppedBy, cut) => {
        if (!final) return end(doneEvent({ text: '', stoppedBy: 'empty', cut, reason: EMPTY_REASONS[why] }));
        const text = stalled ? EMPTY_STALLED : EMPTY_FINAL;
        return end(doneEvent({
            text,
            outcome: {
                status: 'not achieved',
                note: stalled ? 'the analyses stopped adding anything new and no answer was written' : 'no answer was written within the analysis budget',
                body: text
            },
            stoppedBy,
            cut
        }), text);
    };

    const onReply = ({ text, calls, looped, cut }, outOfTime) => {
        if (calls.length) {
            if (!toolsOn || final) invalid('a reply with calls to a request that allowed none (tools off, or the last request)');
            phase = 'running';
            return { next: 'run' };
        }
        const blank = !text.trim();
        const why = looped ? 'loop' : (cut ? 'cut' : 'empty');
        // A reasoning-only, empty or looping reply: one reminder, then ask
        // again (warmer after a loop, and for the rest of the turn).
        if (blank && !reminded) {
            reminded = true;
            warmer = warmer || looped;
            return ask('reminder', outOfTime, { notes: [emptyNudge(why)] });
        }
        const stoppedBy = final ? (stalled ? 'no_progress' : 'max_rounds') : 'answer';
        // An answer is judged as it will read: a conclusion whose analysis
        // did not run may go back to the model while tools can still run in
        // the next round. No answer: the one the gate withdrew stands.
        const judged = blank ? withdrawn : text;
        if (judged === null) return emptyEnding(why, stoppedBy, cut);
        const canSendBack = !blank && toolsOn && !final && !outOfTime();
        const verdict = judge(judged, { canSendBack, rounds: round, stoppedBy, cut });
        if (verdict == null) return emptyEnding(why, stoppedBy, cut);
        if (!isObject(verdict)) invalid('judge must return { sendBack }, { done } or null');
        if (verdict.sendBack) {
            if (!canSendBack) invalid('judge sent an answer back when no send-back was allowed');
            withdrawn = text;
            return ask('gate', outOfTime, { sendBack: verdict.sendBack });
        }
        if (!isObject(verdict.done)) invalid('judge must return { sendBack }, { done } or null');
        // rawText only when the references changed the text (the verification's).
        const { text: answer = judged, rawText, outcome, audit, references, gates } = verdict.done;
        return end(doneEvent({ text: answer, rawText, outcome, stoppedBy, cut, audit, references, gates }));
    };

    const onRound = ({ progressed, waitedMs }, outOfTime) => {
        checkWait(waitedMs);
        waited += waitedMs;
        round += 1;
        noProgressRounds = progressed ? 0 : noProgressRounds + 1;
        if (noProgressRounds >= maxNoProgressRounds) stalled = true;
        return ask('tools', outOfTime);
    };

    const onNoTools = (outOfTime) => {
        if (!toolsOn || round !== 0) invalid('no_tools only answers a request with tools before any round of tools');
        toolsOn = false;
        return ask('no_tools', outOfTime);
    };

    const next = (event) => {
        if (!isObject(event) || typeof event.kind !== 'string') invalid('an event must be { kind, … }');
        if (!EXPECTED[phase].includes(event.kind)) {
            invalid(`a '${event.kind}' event is out of order: ${phase === 'ended' ? 'the turn has ended' : `expected ${EXPECTED[phase].join(' or ')}`}`);
        }
        // The turn's time, read at most once per event (and only when it
        // decides something: never in round 0).
        let elapsed = null;
        const outOfTime = () => {
            if (round === 0) return false;
            if (elapsed === null) elapsed = clock.now() - startedAt - waited;
            return elapsed > maxTurnMs;
        };
        switch (event.kind) {
            case 'start':
                checkWait(event.waitedMs);
                waited += event.waitedMs;
                return ask('first', outOfTime);
            case 'reply': {
                const { text, calls, looped, cut } = event;
                if (typeof text !== 'string' || !Array.isArray(calls) || typeof looped !== 'boolean' || typeof cut !== 'boolean') {
                    invalid('a reply must be { text: string, calls: array, looped: boolean, cut: boolean }');
                }
                return onReply(event, outOfTime);
            }
            case 'round':
                if (typeof event.progressed !== 'boolean') invalid('a round must be { progressed: boolean, waitedMs }');
                return onRound(event, outOfTime);
            default:
                return onNoTools(outOfTime);
        }
    };

    return Object.freeze({ next });
};
