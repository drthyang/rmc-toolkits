// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

import React, { useCallback, useEffect, useMemo, useRef, useState } from 'react';
import { contextToJson } from '../context/runContext';
import { stepSummary } from '../copilot/turnState';
import { useCopilot } from '../useCopilot';
import MarkdownMessage from './MarkdownMessage';
import QuestionCard from './QuestionCard';
import ToolStep from './ToolStep';
import { ArrowUpIcon, ChevronIcon, StopIcon, UserGlyph, WaveMark } from './icons';

// The Copilot conversation: questions, and for each one the assistant's turn —
// what it said, the analyses it ran (ToolStep cards), the clarifying question
// it asked, if any (QuestionCard: the turn waits for the reply), its answer,
// and the Outcome chip with its own verdict. User messages stay plain text; replies
// render Markdown through the sanitizing pipeline (MarkdownMessage), with the
// cards they cite as chips that bring the card on screen (`focus.show`).

const SUGGESTIONS = [
    'Which site is most anharmonic, and does it have a preferred displacement direction?',
    'Compute the bond angles around the central atom of the main polyhedron and compare them with the ideal shape.',
    'Has the modeling run converged, and is the configuration well sampled?'
];

const OUTCOME_TEXT = { achieved: 'Achieved', partial: 'Partial', 'not verified': 'Not verified', 'not achieved': 'Not achieved' };

const STOP_NOTES = {
    stopped: 'Stopped.',
    error: 'Stopped by an error.',
    max_rounds: 'Answered after using the analysis budget for one question.',
    no_progress: 'Answered when repeated analyses stopped adding anything new.',
    empty: 'The model returned an empty reply.'
};

// The stop note, with the reason when the loop gave one (an empty reply).
const stopNoteFor = (turn) => {
    const note = STOP_NOTES[turn.stoppedBy];
    return note && turn.reason ? `${note.replace(/\.$/, '')}: ${turn.reason}.` : note;
};

const formatThinkTime = (ms) => {
    if (!Number.isFinite(ms)) return null;
    const seconds = ms / 1000;
    return seconds < 1 ? '<1s' : `${Math.round(seconds)}s`;
};

// The animated dots shown while the model has started but produced nothing yet.
const ThinkingDots = () => (
    <div className="llm-thinking" role="status" aria-label="Thinking">
        <span className="llm-thinking-dot" />
        <span className="llm-thinking-dot" />
        <span className="llm-thinking-dot" />
    </div>
);

// A collapsible chain-of-thought panel: auto-expanded and shimmering while the
// model is still thinking, then collapsed to "Thought for Ns" (re-expandable).
const ReasoningPanel = ({ reasoning, active, reasoningMs }) => {
    const [manualOpen, setManualOpen] = useState(null);
    const open = manualOpen ?? active;
    const bodyRef = useRef(null);

    useEffect(() => {
        if (active && open && bodyRef.current) bodyRef.current.scrollTop = bodyRef.current.scrollHeight;
    }, [reasoning, active, open]);

    const label = active
        ? 'Thinking'
        : (reasoningMs != null ? `Thought for ${formatThinkTime(reasoningMs)}` : 'Reasoning');

    return (
        <div className={`llm-reasoning${open ? ' is-open' : ''}`}>
            <button
                type="button"
                className="llm-reasoning-toggle"
                onClick={() => setManualOpen(!open)}
                aria-expanded={open}
            >
                <span className={`llm-reasoning-label${active ? ' is-active' : ''}`}>{label}</span>
                <span className="llm-reasoning-chevron" aria-hidden="true"><ChevronIcon /></span>
            </button>
            {open && <div className="llm-reasoning-body" ref={bodyRef}>{reasoning}</div>}
        </div>
    );
};

const TextSegment = ({ item, last, cardLinks }) => {
    const thinkingOnly = item.streaming && !item.text;
    if (!item.text && !item.reasoning && !(item.streaming && last)) return null;
    return (
        <>
            {item.reasoning ? (
                <ReasoningPanel reasoning={item.reasoning} active={thinkingOnly} reasoningMs={item.reasoningMs} />
            ) : null}
            {thinkingOnly && !item.reasoning && <ThinkingDots />}
            {item.text ? (
                <div className={`llm-bubble is-markdown${item.streaming ? ' is-streaming' : ''}`}>
                    <MarkdownMessage {...cardLinks}>{item.text}</MarkdownMessage>
                    {item.streaming && <span className="llm-caret" aria-hidden="true" />}
                </div>
            ) : null}
        </>
    );
};

// Consecutive tool steps render as one list between text segments.
const groupItems = (items) => items.reduce((groups, item) => {
    const previous = groups[groups.length - 1];
    if (item.kind === 'tool' && previous?.kind === 'steps') previous.steps.push(item);
    else if (item.kind === 'tool') groups.push({ kind: 'steps', steps: [item] });
    else groups.push(item);
    return groups;
}, []);

const AssistantTurn = ({ turn, pages, onOpenPage, cardLinks, question = null, onAnswer = null }) => {
    const groups = groupItems(turn.items);
    const summary = stepSummary(turn.items);
    const stopNote = stopNoteFor(turn);
    // An empty reply's reason already says the length limit was reached.
    const cutNote = turn.cut && !(turn.stoppedBy === 'empty' && turn.reason);
    return (
        <div className="llm-message is-assistant">
            <span className="llm-message-avatar" aria-hidden="true"><WaveMark /></span>
            <div className="llm-message-body llm-turn">
                <span className="llm-message-name">Copilot</span>
                {groups.map((group, index) => (group.kind === 'steps' ? (
                    <ol key={index} className="llm-steps" aria-label="Analyses run">
                        {group.steps.map((step) => (step.interactive ? (
                            <QuestionCard
                                key={step.id}
                                step={step}
                                request={step.status === 'running' ? question : null}
                                onAnswer={onAnswer}
                            />
                        ) : (
                            <ToolStep
                                key={step.id}
                                step={step}
                                pageLabel={step.page ? pages?.[step.page] : null}
                                onOpenPage={onOpenPage}
                            />
                        )))}
                    </ol>
                ) : (
                    <TextSegment key={index} item={group} last={turn.streaming && index === groups.length - 1} cardLinks={cardLinks} />
                )))}
                {turn.audit?.length > 0 && (
                    // Checks on the answer itself (answerAudit.js): where it
                    // contradicts what the tools or the run context established.
                    <ul className="llm-audit" aria-label="Checks on this answer">
                        {turn.audit.map((finding) => (
                            <li key={finding.rule} className="llm-audit-item">
                                <span className="llm-audit-label">Check</span>
                                {finding.text}
                            </li>
                        ))}
                    </ul>
                )}
                {turn.outcome && (
                    <div className={`llm-outcome is-${turn.outcome.status.replace(' ', '-')}`} role="status">
                        <span className="llm-outcome-label">Outcome</span>
                        <strong>{OUTCOME_TEXT[turn.outcome.status]}</strong>
                        {turn.outcome.note && (
                            // Models write math and card links in the note too.
                            <span className="llm-outcome-note"><MarkdownMessage inline {...cardLinks}>{turn.outcome.note}</MarkdownMessage></span>
                        )}
                    </div>
                )}
                {!turn.streaming && (summary || stopNote || turn.toolsUnsupported || turn.cut || turn.gates?.length > 0) && (
                    <p className={`llm-turn-footer${turn.cut ? ' is-cut' : ''}`}>
                        {[
                            summary,
                            // The evidence gate's send-backs (THY-18).
                            ...(turn.gates || []).map((gate) => gate.notice),
                            turn.toolsUnsupported ? 'this model cannot run analyses — answered from the run summary' : '',
                            cutNote ? 'the reply reached the model\'s length limit and is cut short' : '',
                            turn.stoppedBy !== 'answer' ? stopNote : ''
                        ].filter(Boolean).join(' · ')}
                    </p>
                )}
            </div>
        </div>
    );
};

// One short line for screen readers when a step starts or ends, or a turn
// ends with its verdict — never the streamed tokens themselves.
const announcementFor = (turn) => {
    // The newest turn can be the question itself (a turn that failed before
    // any output is not kept): nothing to announce.
    if (!turn || !Array.isArray(turn.items)) return '';
    if (!turn.streaming && turn.outcome) {
        const checks = (turn.audit || []).map((finding) => ` Check: ${finding.text}`).join('');
        return `Outcome: ${OUTCOME_TEXT[turn.outcome.status]}. ${turn.outcome.note || ''}${checks}`.trim();
    }
    const step = [...turn.items].reverse().find((item) => item.kind === 'tool');
    if (!step) return turn.streaming ? '' : 'Answer ready.';
    if (step.interactive) {
        const asked = step.args?.question || '';
        return step.status === 'running' ? `The Copilot asks: ${asked}` : `Answered: ${step.data?.answer ?? 'skipped'}`;
    }
    return step.status === 'running' ? `Running: ${step.label}` : `${step.label}: ${step.status}`;
};

const UserMessage = ({ content }) => (
    <div className="llm-message is-user">
        <span className="llm-message-avatar" aria-hidden="true"><UserGlyph /></span>
        <div className="llm-message-body">
            <span className="llm-message-name">You</span>
            <div className="llm-bubble">{content}</div>
        </div>
    </div>
);

const CopilotView = ({
    context, run, settings, connected, tools, toolsEnabled, pages, onOpenPage, runKey,
    cards = null, focus = null, onTurnEnd = null, onBusyChange = null, onQuestion = null, skills = null
}) => {
    const copilot = useCopilot({ settings, run, tools, toolsEnabled, runKey, focus, onTurnEnd, onQuestion, skills });
    useEffect(() => { onBusyChange?.(copilot.busy); }, [copilot.busy, onBusyChange]);
    const showCard = useCallback((id) => { focus?.show([id], { source: 'user' }); }, [focus]);
    const previewCard = useCallback((id) => { focus?.preview(id); }, [focus]);
    const cardLinks = useMemo(() => (cards?.length && focus
        ? { cards, onCardClick: showCard, onCardHover: previewCard }
        : null), [cards, focus, showCard, previewCard]);
    const [draft, setDraft] = useState('');
    const logRef = useRef(null);
    const textareaRef = useRef(null);
    const hasConversation = copilot.turns.length > 0 || copilot.busy;
    const announcement = useMemo(
        () => announcementFor(copilot.live || copilot.turns.at(-1)),
        [copilot.live, copilot.turns]
    );

    // Keep the newest message in view as turns, steps and tokens arrive — but
    // never scroll a waiting question's top out of view (a short sheet on a
    // phone would show its options without the question).
    useEffect(() => {
        const el = logRef.current;
        if (!el) return;
        const waiting = el.querySelector('.llm-ask.is-waiting');
        if (!waiting) {
            el.scrollTop = el.scrollHeight;
            return;
        }
        const top = waiting.getBoundingClientRect().top - el.getBoundingClientRect().top + el.scrollTop;
        el.scrollTop = Math.max(0, Math.min(el.scrollHeight - el.clientHeight, top - 8));
    }, [copilot.turns, copilot.live, copilot.question]);

    // The chatbox focuses its panel when it opens before the model is
    // connected (the composer is disabled until then); once it connects,
    // move focus on to the composer — unless the user has moved it elsewhere.
    useEffect(() => {
        const ta = textareaRef.current;
        if (!connected || !ta) return;
        const panel = ta.closest('.llm-chat-panel');
        if (panel && document.activeElement === panel) ta.focus();
    }, [connected]);

    // A turn that ends with the focus lost (the question card it was on
    // re-rendered as the finished turn) hands it to the composer, not to the
    // page's start.
    const wasBusy = useRef(false);
    useEffect(() => {
        if (copilot.busy) {
            wasBusy.current = true;
            return;
        }
        if (!wasBusy.current) return;
        wasBusy.current = false;
        const ta = textareaRef.current;
        const active = document.activeElement;
        if (ta && !ta.disabled && (!active || active === document.body)) ta.focus({ preventScroll: true });
    }, [copilot.busy]);

    const resizeTextarea = () => {
        const ta = textareaRef.current;
        if (!ta) return;
        ta.style.height = 'auto';
        ta.style.height = `${Math.min(ta.scrollHeight, 160)}px`;
    };

    const send = () => {
        const question = draft.trim();
        if (!question || copilot.busy || !connected) return;
        setDraft('');
        requestAnimationFrame(resizeTextarea);
        copilot.ask(question);
    };

    const applySuggestion = (text) => {
        setDraft(text);
        const ta = textareaRef.current;
        if (ta) {
            ta.focus();
            requestAnimationFrame(resizeTextarea);
        }
    };

    // What the Copilot can run on this run with its analyses on (the turn's
    // own ask_user and read_skill included): the status line and the list.
    const offered = copilot.offered;
    const toolStatus = (() => {
        if (copilot.toolsUnsupported) return 'This model cannot run analyses — answers come from the run summary.';
        if (!toolsEnabled) return 'Analyses are off — answers come from the run summary.';
        const analyses = offered.tools.filter((tool) => !tool.view && !tool.interactive && !tool.keep).length;
        if (!analyses) return 'No analyses are available for this run.';
        const pointing = offered.view ? ' and point at its cards' : '';
        const count = offered.skills.length;
        const methods = count ? `, following ${count} skill${count === 1 ? '' : 's'}` : '';
        return `Can run ${analyses} ${analyses === 1 ? 'analysis' : 'analyses'} on this run${pointing}${methods}.`;
    })();

    return (
        <div className="llm-chat">
            {hasConversation ? (
                <div className="llm-chat-log" ref={logRef}>
                    {copilot.turns.map((turn, index) => (turn.role === 'user'
                        ? <UserMessage key={index} content={turn.content} />
                        : <AssistantTurn key={index} turn={turn} pages={pages} onOpenPage={onOpenPage} cardLinks={cardLinks} />))}
                    {copilot.live && (
                        <AssistantTurn
                            turn={copilot.live}
                            pages={pages}
                            onOpenPage={onOpenPage}
                            cardLinks={cardLinks}
                            question={copilot.question}
                            onAnswer={copilot.answerQuestion}
                        />
                    )}
                </div>
            ) : (
                <div className="llm-chat-empty">
                    <span className="llm-chat-empty-icon" aria-hidden="true"><WaveMark /></span>
                    <p className="llm-chat-empty-title">What should I check?</p>
                    <p className="llm-chat-empty-sub">
                        The Copilot runs the Workbench analyses itself, checks what comes back, and tells you whether
                        the answer holds.
                    </p>
                    {connected && (
                        <div className="llm-suggestions">
                            {SUGGESTIONS.map((text) => (
                                <button key={text} type="button" className="llm-suggestion" onClick={() => applySuggestion(text)}>
                                    {text}
                                </button>
                            ))}
                        </div>
                    )}
                    {offered.tools.length > 0 && (
                        <details className="llm-context-inspector">
                            <summary title="The analyses the Copilot may run on this run, and its view tools. Analyses only compute; view tools only point at cards and turn 3D views — none changes data or a file.">
                                What the Copilot can run
                            </summary>
                            <ul className="llm-tool-list">
                                {offered.tools.map((tool) => (
                                    <li key={tool.name}><code>{tool.name}</code> — {tool.description}</li>
                                ))}
                            </ul>
                        </details>
                    )}
                    {context && (
                        <details className="llm-context-inspector">
                            <summary title="The run's metrics, symmetry and convergence history travel with every question, so answers can quote the actual numbers.">
                                Context sent to the model
                            </summary>
                            <pre>{contextToJson(context)}</pre>
                        </details>
                    )}
                </div>
            )}

            <div className="llm-sr-only" aria-live="polite">{announcement}</div>
            {copilot.error && <div className="llm-connection-status is-error" role="alert">{copilot.error}</div>}
            {copilot.skillsError && <div className="llm-connection-status is-error" role="alert">The Copilot's skills could not be read: {copilot.skillsError}</div>}

            <form
                className="llm-composer"
                onSubmit={(event) => {
                    event.preventDefault();
                    send();
                }}
            >
                <textarea
                    ref={textareaRef}
                    className="llm-composer-input"
                    rows={1}
                    value={draft}
                    placeholder={connected ? 'Ask the Copilot to check something in this run…' : 'Connect a model to start…'}
                    aria-label="Message the Copilot"
                    onChange={(event) => {
                        setDraft(event.target.value);
                        resizeTextarea();
                    }}
                    onKeyDown={(event) => {
                        if (event.key === 'Enter' && !event.shiftKey) {
                            event.preventDefault();
                            send();
                        }
                    }}
                    disabled={!connected}
                />
                {copilot.busy ? (
                    <button type="button" className="llm-send is-stop" onClick={copilot.stop} aria-label="Stop">
                        <StopIcon />
                    </button>
                ) : (
                    <button type="submit" className="llm-send" disabled={!connected || !draft.trim()} aria-label="Send message">
                        <ArrowUpIcon />
                    </button>
                )}
            </form>
            <div className="llm-composer-foot">
                <span className="llm-composer-status">{toolStatus}</span>
                {hasConversation && !copilot.busy && (
                    <button type="button" className="llm-link-btn" onClick={copilot.clear}>New conversation</button>
                )}
            </div>
        </div>
    );
};

export default CopilotView;
