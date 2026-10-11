// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

import React, { useEffect, useId, useRef, useState } from 'react';

// The Copilot's clarifying question (ask_user, copilot/askUser.js) as a card
// in the conversation: the question, 2–4 option buttons (each with its one-line
// meaning; the Copilot's default is marked), an "Other…" field and Skip. The
// turn waits until one of them answers. Keyboard: Tab through the options,
// ↑/↓ (or ←/→) move between them, 1–4 pick one, Enter in the field sends it.
// When the card appears it takes focus, unless the user is typing elsewhere.
// Once answered (or stopped) it keeps the question and the reply, compact.

// Focus may move to the card from nowhere in particular or from inside the
// conversation (the empty composer, where it rests after a question is sent),
// never out of a field that holds text the user is typing.
const mayTakeFocus = (card) => {
    if (typeof document === 'undefined' || !card) return false;
    const active = document.activeElement;
    if (!active || active === document.body) return true;
    if (active.matches?.('input, textarea, select, [contenteditable="true"]') && (active.value || active.isContentEditable)) {
        return false;
    }
    return Boolean(card.closest('.llm-chat-panel')?.contains(active));
};

const repliedText = (step) => {
    const data = step.data || {};
    if (data.skipped) return data.answer ? `Skipped: the Copilot uses its default, ${data.answer}.` : 'Skipped.';
    return data.answer ?? '';
};

const QuestionCard = ({ step, request = null, onAnswer }) => {
    const questionId = useId();
    const otherId = useId();
    const cardRef = useRef(null);
    const optionRefs = useRef([]);
    const [other, setOther] = useState('');
    const waiting = step.status === 'running' && Boolean(request);
    const question = request?.question || step.data?.question || step.args?.question || '';

    useEffect(() => {
        if (!waiting) return;
        const target = optionRefs.current[request.defaultIndex ?? 0] || optionRefs.current[0];
        if (target && mayTakeFocus(cardRef.current)) target.focus({ preventScroll: true });
    }, [waiting, request]);

    const choose = (index) => onAnswer?.({ answer: request.options[index].label, source: 'option', index });
    const sendOther = () => {
        const text = other.trim();
        if (text) onAnswer?.({ answer: text, source: 'other', index: null });
    };

    const onOptionKey = (event, index) => {
        const count = request.options.length;
        const move = { ArrowDown: 1, ArrowRight: 1, ArrowUp: -1, ArrowLeft: -1 }[event.key];
        if (move) {
            event.preventDefault();
            optionRefs.current[(index + move + count) % count]?.focus();
            return;
        }
        const digit = Number(event.key);
        if (Number.isInteger(digit) && digit >= 1 && digit <= count) {
            event.preventDefault();
            choose(digit - 1);
        }
    };

    if (!waiting) {
        const stopped = step.status === 'stopped' || (step.status === 'running' && !request);
        const failed = step.status === 'failed';
        return (
            <li className={`llm-step llm-ask is-${stopped ? 'stopped' : failed ? 'failed' : 'answered'}`} ref={cardRef}>
                <div className="llm-ask-kicker">
                    <span className="llm-step-dot" aria-hidden="true" />
                    <span>Question for you</span>
                    <span className="llm-step-meta">
                        {step.status === 'running' ? 'Waiting…' : stopped ? 'Not answered' : failed ? 'Not asked' : 'Answered'}
                    </span>
                </div>
                {question && <p className="llm-ask-question">{question}</p>}
                {failed && step.error && <p className="llm-step-headline is-fail">{step.error}</p>}
                {!stopped && !failed && step.status !== 'running' && (
                    <p className="llm-ask-reply"><span className="llm-ask-reply-who">You</span>{repliedText(step)}</p>
                )}
            </li>
        );
    }

    return (
        <li
            className="llm-step llm-ask is-waiting"
            ref={cardRef}
            role="group"
            aria-labelledby={questionId}
        >
            <div className="llm-ask-kicker">
                <span className="llm-step-dot" aria-hidden="true" />
                <span>Question for you</span>
                <span className="llm-step-meta">The Copilot waits for your answer</span>
            </div>
            <p className="llm-ask-question" id={questionId}>{question}</p>
            <div className="llm-ask-options" role="group" aria-label="Choices">
                {request.options.map((option, index) => (
                    <button
                        key={option.label}
                        ref={(element) => { optionRefs.current[index] = element; }}
                        type="button"
                        className={`llm-ask-option${index === request.defaultIndex ? ' is-default' : ''}`}
                        onClick={() => choose(index)}
                        onKeyDown={(event) => onOptionKey(event, index)}
                        aria-keyshortcuts={String(index + 1)}
                    >
                        <span className="llm-ask-option-key" aria-hidden="true">{index + 1}</span>
                        <span className="llm-ask-option-text">
                            <span className="llm-ask-option-label">
                                {option.label}
                                {index === request.defaultIndex && <span className="llm-ask-default">default</span>}
                            </span>
                            {option.description && <span className="llm-ask-option-desc">{option.description}</span>}
                        </span>
                    </button>
                ))}
            </div>
            {request.allowOther && (
                <form
                    className="llm-ask-other"
                    onSubmit={(event) => {
                        event.preventDefault();
                        sendOther();
                    }}
                >
                    <label className="llm-sr-only" htmlFor={otherId}>Another answer</label>
                    <input
                        id={otherId}
                        className="llm-input llm-ask-other-input"
                        type="text"
                        value={other}
                        placeholder="Other…"
                        autoComplete="off"
                        onChange={(event) => setOther(event.target.value)}
                    />
                    <button type="submit" className="llm-primary llm-ask-send" disabled={!other.trim()}>Send</button>
                </form>
            )}
            <button
                type="button"
                className="llm-link-btn llm-ask-skip"
                onClick={() => onAnswer?.({ answer: null, source: 'skipped', index: null })}
                title={request.defaultIndex != null
                    ? 'Let the Copilot go on with its default'
                    : 'Let the Copilot choose; it says what it assumed'}
            >
                Skip
            </button>
        </li>
    );
};

export default QuestionCard;
