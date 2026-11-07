// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

import React, { useCallback, useEffect, useId, useMemo, useRef, useState } from 'react';
import { contextToJson } from '../context/runContext';
import { effectOf } from '../copilot/effects';
import { STOPPED_BY } from '../copilot/turnPolicy';
import { changeReply, questionReply, stepSummary, turnReferences } from '../copilot/turnState';
import { useCopilot } from '../useCopilot';
import ChangeCard from './ChangeCard';
import ContextMeter from './ContextMeter';
import { compactedNote, memoryNote } from './chatNotes';
import MarkdownMessage from './MarkdownMessage';
import { answerParts } from './answerParts';
import { plainReferences } from './referenceDisplay';
import QuestionCard from './QuestionCard';
import ToolStep from './ToolStep';
import { ArrowUpIcon, ChevronIcon, StopIcon, UserGlyph, WaveMark } from './icons';

// The Copilot conversation: questions, and for each one the assistant's turn —
// what it said, the analyses it ran (ToolStep cards), the clarifying question
// it asked, if any (QuestionCard: the turn waits for the reply), a change it
// proposed to a page's inputs (ChangeCard: Apply / Decline, then Undo), its
// answer, and the Outcome chip with its own verdict. User messages stay plain text; replies
// render Markdown through the sanitizing pipeline (MarkdownMessage), with the
// cards they cite as chips that bring the card on screen (`focus.show`).

const SUGGESTIONS = [
    'Which site is most anharmonic, and does it have a preferred displacement direction?',
    'Compute the bond angles around the central atom of the main polyhedron and compare them with the ideal shape.',
    'Has the modeling run converged, and is the configuration well sampled?'
];

const OUTCOME_TEXT = { achieved: 'Achieved', partial: 'Partial', 'not verified': 'Not verified', 'not achieved': 'Not achieved' };

// The turn's ending as its note (turnPolicy.js STOPPED_BY: none for an
// answer), with the reason when the turn gave one (an empty reply).
const stopNoteFor = (turn) => {
    const note = STOPPED_BY[turn.stoppedBy]?.note;
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

// `references`: the turn's (turnReferences), so a reference the turn's
// results fill shows its number while the text still streams. A long
// finished text shows its lead — the conclusion, which the prompt asks for
// first — with the rest behind "Show details" (answerParts.js); a text
// still streaming shows whole.
const TextSegment = ({ item, last, cardLinks, references }) => {
    const [open, setOpen] = useState(false);
    const detailsId = useId();
    const parts = useMemo(() => (item.streaming || !item.text ? null : answerParts(item.text)), [item.streaming, item.text]);
    const thinkingOnly = item.streaming && !item.text;
    if (!item.text && !item.reasoning && !(item.streaming && last)) return null;
    const markdown = (text) => (
        <MarkdownMessage {...cardLinks} references={references} streaming={Boolean(item.streaming)}>{text}</MarkdownMessage>
    );
    return (
        <>
            {item.reasoning ? (
                <ReasoningPanel reasoning={item.reasoning} active={thinkingOnly} reasoningMs={item.reasoningMs} />
            ) : null}
            {thinkingOnly && !item.reasoning && <ThinkingDots />}
            {item.text ? (
                <div className={`llm-bubble is-markdown${item.streaming ? ' is-streaming' : ''}${parts?.rest ? ' is-folded' : ''}`}>
                    {parts?.rest ? (
                        <>
                            {markdown(parts.lead)}
                            <button
                                type="button"
                                className="llm-answer-toggle"
                                aria-expanded={open}
                                aria-controls={detailsId}
                                onClick={() => setOpen((value) => !value)}
                            >
                                {open ? 'Hide details' : 'Show details'}
                                <span className="llm-answer-chevron" aria-hidden="true"><ChevronIcon /></span>
                            </button>
                            <div id={detailsId} className="llm-answer-details" hidden={!open}>
                                {open ? markdown(parts.rest) : null}
                            </div>
                        </>
                    ) : markdown(item.text)}
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

const SUMMARY_DIVIDER_TEXT = 'Summarized — the answers above go to the model as their sentences with numbers';

const SummaryDivider = () => (
    <p className="llm-config-divider" title="Summarize and continue: each earlier answer is re-sent as its sentences with numbers and its Outcome line, at most 600 characters. No model writes the summary.">
        <span>{SUMMARY_DIVIDER_TEXT}</span>
    </p>
);

// One assistant turn as the dock draws it (exported for the tests that render it alone).
export const AssistantTurn = ({
    turn, pages, onOpenPage, cardLinks, question = null, onAnswer = null, approval = null, onDecide = null, undos = null,
    onUndo = null
}) => {
    const groups = groupItems(turn.items);
    const summary = stepSummary(turn.items);
    // Filled numbers while streaming; never a raw {{…}} (referenceDisplay.js).
    const references = useMemo(() => turnReferences(turn.items, turn.shown), [turn.items, turn.shown]);
    const stopNote = stopNoteFor(turn);
    // An empty reply's reason already says the length limit was reached.
    const cutNote = turn.cut && !(turn.stoppedBy === 'empty' && turn.reason);
    const stepCard = (step) => {
        if (step.interactive) {
            return (
                <QuestionCard
                    key={step.id}
                    step={step}
                    request={step.status === 'running' ? question : null}
                    onAnswer={onAnswer}
                />
            );
        }
        if (step.effect === 'change') {
            return (
                <ChangeCard
                    key={step.id}
                    step={step}
                    request={step.status === 'running' ? approval : null}
                    onDecide={onDecide}
                    undo={(step.change?.id && undos?.[step.change.id]) || null}
                    onUndo={onUndo}
                    pageLabel={step.page ? pages?.[step.page] : null}
                    onOpenPage={onOpenPage}
                />
            );
        }
        return (
            <ToolStep
                key={step.id}
                step={step}
                pageLabel={step.page ? pages?.[step.page] : null}
                onOpenPage={onOpenPage}
            />
        );
    };
    return (
        <div className="llm-message is-assistant">
            <span className="llm-message-avatar" aria-hidden="true"><WaveMark /></span>
            <div className="llm-message-body llm-turn">
                <span className="llm-message-name">Copilot</span>
                {turn.memory && <p className="llm-memory-note">{memoryNote(turn.memory)}</p>}
                {groups.map((group, index) => (group.kind === 'steps' ? (
                    <ol key={index} className="llm-steps" aria-label="Analyses run">
                        {group.steps.map(stepCard)}
                    </ol>
                ) : (
                    <TextSegment key={index} item={group} last={turn.streaming && index === groups.length - 1} cardLinks={cardLinks} references={references} />
                )))}
                {/* The Outcome, then the checks its note points to ("a check below"). */}
                {turn.outcome && (
                    <div className={`llm-outcome is-${turn.outcome.status.replace(' ', '-')}`} role="status">
                        <span className="llm-outcome-label">Outcome</span>
                        <strong>{OUTCOME_TEXT[turn.outcome.status]}</strong>
                        {turn.outcome.note && (
                            // Models write math and card links in the note too.
                            <span className="llm-outcome-note"><MarkdownMessage inline {...cardLinks} references={references}>{turn.outcome.note}</MarkdownMessage></span>
                        )}
                    </div>
                )}
                {turn.audit?.length > 0 && (
                    // Checks on the answer itself (the verification, facts.js): where it
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
                {!turn.streaming && (summary || stopNote || turn.toolsUnsupported || turn.cut || turn.gates?.length > 0 || turn.compacted > 0) && (
                    <p className={`llm-turn-footer${turn.cut ? ' is-cut' : ''}`}>
                        {[
                            summary,
                            // The evidence gate's send-backs (THY-18).
                            ...(turn.gates || []).map((gate) => gate.notice),
                            turn.toolsUnsupported ? 'this model cannot run analyses — answered from the run summary' : '',
                            // Never silently: earlier results sent shortened (agentLoop.js turnCharBudget).
                            compactedNote(turn.compacted),
                            cutNote ? 'the reply reached the model\'s length limit and is cut short' : '',
                            stopNote
                        ].filter(Boolean).join(' · ')}
                    </p>
                )}
            </div>
        </div>
    );
};

// A question step's line, in the card's words (questionReply): asked, then
// answered, skipped (the default stands), not asked, or not answered.
const QUESTION_LINE = {
    answered: (text) => `Answered: ${text}`,
    skipped: (text) => text,
    not_asked: (text) => `Not asked: ${text}`,
    stopped: () => 'Not answered.'
};

// A change step's line, in the card's words (changeReply).
const CHANGE_LINE = {
    applied: (text) => `Change applied: ${text}`,
    declined: (text) => `Change declined: ${text}`,
    unchanged: (text) => `No change needed: ${text}`,
    stale: () => 'Change not applied: the page changed while it waited.',
    not_applied: (text) => `Change not applied: ${text}`,
    stopped: () => 'Change stopped before it was made.'
};

// One short line for screen readers when a step starts or ends, or a turn
// ends with its verdict — never the streamed tokens themselves.
const announcementFor = (turn, approval = null) => {
    // The newest turn can be the question itself (a turn that failed before
    // any output is not kept): nothing to announce.
    if (!turn || !Array.isArray(turn.items)) return '';
    if (!turn.streaming && turn.outcome) {
        const checks = (turn.audit || []).map((finding) => ` Check: ${finding.text}`).join('');
        return `Outcome: ${OUTCOME_TEXT[turn.outcome.status]}. ${plainReferences(turn.outcome.note || '')}${checks}`.trim();
    }
    const step = [...turn.items].reverse().find((item) => item.kind === 'tool');
    if (!step) return turn.streaming ? '' : 'Answer ready.';
    if (step.interactive) {
        if (step.status === 'running') return `The Copilot asks: ${step.args?.question || ''}`;
        const { state, text } = questionReply(step);
        return QUESTION_LINE[state](text);
    }
    if (step.effect === 'change') {
        if (step.status === 'running') return approval ? `The Copilot asks to apply a change: ${approval.line}` : `Running: ${step.label}`;
        const { state, text } = changeReply(step);
        return CHANGE_LINE[state]?.(text) ?? `${step.label}: ${step.status}`;
    }
    return step.status === 'running' ? `Running: ${step.label}` : `${step.label}: ${step.status}`;
};

// The divider the dock draws where a new configuration of the run was saved
// (configKey, useCopilot.js configurationBreaks): the answers above it are
// not sent to the model with the next question.
const CONFIGURATION_DIVIDER_TEXT = 'New configuration saved — earlier answers were about the previous one';

const ConfigurationDivider = () => (
    <p
        className="llm-config-divider"
        title="The run's .rmc6f changed. The answers above are kept on screen but are not sent to the model with the next question."
    >
        <span>{CONFIGURATION_DIVIDER_TEXT}</span>
    </p>
);

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
    context, run, settings, connected, tools, toolsEnabled, pages, onOpenPage, runKey, configKey = null,
    cards = null, focus = null, onTurnEnd = null, onBusyChange = null, onQuestion = null, onApproval = null, skills = null,
    autoApprove = false
}) => {
    const copilot = useCopilot({
        settings, run, tools, toolsEnabled, runKey, configKey, focus, onTurnEnd, onQuestion, onApproval, skills, autoApprove, connected
    });
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
    // The newest turn's line, then the divider's once a new configuration
    // ends the conversation (after the running answer, if any).
    const configurationChanged = copilot.breaks.after;
    const announcement = useMemo(() => {
        const line = announcementFor(copilot.live || copilot.turns.at(-1), copilot.approval);
        return configurationChanged && !copilot.live ? `${line} ${CONFIGURATION_DIVIDER_TEXT}.`.trim() : line;
    }, [copilot.live, copilot.turns, copilot.approval, configurationChanged]);

    // Keep the newest message in view as turns, steps and tokens arrive — but
    // never scroll a waiting question's or change's top out of view (a short
    // sheet on a phone would show its buttons without what they answer).
    useEffect(() => {
        const el = logRef.current;
        if (!el) return;
        const waiting = el.querySelector('.llm-ask.is-waiting, .llm-change.is-waiting');
        if (!waiting) {
            el.scrollTop = el.scrollHeight;
            return;
        }
        const top = waiting.getBoundingClientRect().top - el.getBoundingClientRect().top + el.scrollTop;
        el.scrollTop = Math.max(0, Math.min(el.scrollHeight - el.clientHeight, top - 8));
    }, [copilot.turns, copilot.live, copilot.question, copilot.approval]);

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
    // Analyses are the compute tools (copilot/effects.js) that are neither a
    // question nor a skill read; a change tool sets a page's inputs once the
    // user approves (Auto: without asking).
    const offered = copilot.offered;
    const toolStatus = (() => {
        if (copilot.toolsUnsupported) return 'This model cannot run analyses — answers come from the run summary.';
        if (!toolsEnabled) return 'Analyses are off — answers come from the run summary.';
        const analyses = offered.tools.filter((tool) => effectOf(tool) === 'compute' && !tool.interactive && !tool.keep).length;
        if (!analyses) return 'No analyses are available for this run.';
        const pointing = offered.view ? ' and point at its cards' : '';
        const changes = offered.tools.some((tool) => effectOf(tool) === 'change')
            ? (autoApprove ? '; changes its inputs without asking (Auto)' : '; changes its inputs once you approve')
            : '';
        const count = offered.skills.length;
        const methods = count ? `, following ${count} skill${count === 1 ? '' : 's'}` : '';
        return `Can run ${analyses} ${analyses === 1 ? 'analysis' : 'analyses'} on this run${pointing}${methods}${changes}.`;
    })();

    return (
        <div className="llm-chat">
            {hasConversation ? (
                <div className="llm-chat-log" ref={logRef}>
                    {copilot.turns.flatMap((turn, index) => [
                        ...(copilot.breaks.before.includes(index) ? [<ConfigurationDivider key={`configuration-${index}`} />] : []),
                        ...(copilot.summarizedAt === index ? [<SummaryDivider key={`summary-${index}`} />] : []),
                        turn.role === 'user'
                            ? <UserMessage key={index} content={turn.content} />
                            : (
                                <AssistantTurn
                                    key={index}
                                    turn={turn}
                                    pages={pages}
                                    onOpenPage={onOpenPage}
                                    cardLinks={cardLinks}
                                    undos={copilot.undos}
                                    onUndo={copilot.undoChange}
                                />
                            )
                    ])}
                    {copilot.live && (
                        <AssistantTurn
                            turn={copilot.live}
                            pages={pages}
                            onOpenPage={onOpenPage}
                            cardLinks={cardLinks}
                            question={copilot.question}
                            onAnswer={copilot.answerQuestion}
                            approval={copilot.approval}
                            onDecide={copilot.decideChange}
                            undos={copilot.undos}
                            onUndo={copilot.undoChange}
                        />
                    )}
                    {configurationChanged && <ConfigurationDivider />}
                    {copilot.summarizedAt === copilot.turns.length && !copilot.live && <SummaryDivider />}
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
                            <summary title="The analyses the Copilot may run on this run, its view tools and its change tools. Analyses only compute; view tools only point at cards and turn 3D views; a change tool sets a page's inputs once you approve — none changes data or a file.">
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
                <span className="llm-composer-tail">
                    {connected && (
                        <ContextMeter
                            meter={copilot.meter}
                            busy={copilot.busy}
                            onStartFresh={hasConversation ? copilot.clear : null}
                            onSummarize={copilot.summarize}
                        />
                    )}
                    {hasConversation && !copilot.busy && !(connected && ['amber', 'red'].includes(copilot.meter?.level)) && (
                        <button type="button" className="llm-link-btn" onClick={copilot.clear}>New conversation</button>
                    )}
                </span>
            </div>
        </div>
    );
};

export default CopilotView;
