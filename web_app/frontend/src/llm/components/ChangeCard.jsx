// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

import React, { useEffect, useId, useRef, useState } from 'react';
import { changeReply } from '../copilot/turnState';
import { mayTakeFocus } from './cardFocus';

// A change the Copilot makes to a page's inputs (a change tool,
// copilot/effects.js) as a card in the conversation. While it waits for the
// user (Ask first): one plain line — the page after the change — with each
// control's before → after behind "?", and Apply / Decline; the turn waits.
// Once decided it keeps the line and how it went: applied (the Copilot's
// change, with Undo), declined, not applied (the page changed while the card
// waited, or the change was refused, with why), or already so. In Auto a
// change runs without waiting and its card says so. The state is read from
// the step's result and its change record (turnState.js changeReply), and
// the Undo's from the dock (`undo`: { state: 'running' | 'done' | 'failed',
// reason }). When the card starts waiting it takes the focus, unless the
// user is typing elsewhere (cardFocus.js), so Tab reaches Apply.

const META = {
    running: 'Working…',
    applied: 'Applied by the Copilot',
    declined: 'Declined',
    unchanged: 'No change',
    stale: 'Not applied',
    not_applied: 'Not applied',
    stopped: 'Stopped'
};

// One short line for the user, where the step's own text is the model's.
const SAY = {
    stale: 'Not applied: the page changed while this waited for you.',
    unchanged: 'Already on the page.',
    declined: 'Declined: nothing changed.',
    stopped: 'Stopped before the change was made.'
};

const capitalize = (text) => (text ? `${text.charAt(0).toUpperCase()}${text.slice(1)}` : text);

const UNDO_TEXT = {
    running: 'Undoing…',
    done: 'Undone: the earlier inputs are back.'
};

const ChangeCard = ({ step, request = null, onDecide, undo = null, onUndo, pageLabel = null, onOpenPage }) => {
    const titleId = useId();
    const detailsId = useId();
    const cardRef = useRef(null);
    const [open, setOpen] = useState(false);
    const waiting = step.status === 'running' && Boolean(request);
    const reply = changeReply(step);
    const record = step.change;
    const line = waiting ? request.line : (record?.line || step.label || step.name);
    const details = waiting ? request.details || [] : record?.details || [];
    const kicker = `${pageLabel || 'Page'} change`;

    useEffect(() => {
        if (waiting && mayTakeFocus(cardRef.current)) cardRef.current?.focus({ preventScroll: true });
    }, [waiting]);
    // Decided: the focused button is gone, so the card keeps the focus.
    const wasWaiting = useRef(false);
    useEffect(() => {
        if (waiting) {
            wasWaiting.current = true;
            return;
        }
        if (!wasWaiting.current) return;
        wasWaiting.current = false;
        const active = typeof document !== 'undefined' ? document.activeElement : null;
        if (!active || active === document.body || cardRef.current?.contains(active)) cardRef.current?.focus({ preventScroll: true });
    }, [waiting]);

    const state = waiting ? 'waiting' : reply.state;
    const applied = state === 'applied';
    const undone = applied && undo?.state === 'done';
    const meta = waiting
        ? 'Needs your approval'
        : (undone ? 'Undone' : `${META[state] || state}${applied && record?.auto ? ' · Auto' : ''}`);
    // The tool's short line for the user (the page did not compute, no angles).
    const note = applied ? record?.note : null;
    const said = state === 'not_applied' ? capitalize(reply.text) : SAY[state];

    return (
        <li
            className={`llm-step llm-change is-${state}`}
            ref={cardRef}
            tabIndex={-1}
            role="group"
            aria-labelledby={titleId}
            data-change-state={state}
        >
            <div className="llm-ask-kicker">
                <span className="llm-step-dot" aria-hidden="true" />
                <span>{kicker}</span>
                <span className="llm-step-meta">{meta}</span>
            </div>
            <div className="llm-change-line">
                <p className="llm-change-text" id={titleId} title={waiting ? request.title : undefined}>{line}</p>
                {details.length > 0 && (
                    <button
                        type="button"
                        className="llm-change-help"
                        aria-expanded={open}
                        aria-controls={detailsId}
                        aria-label={open ? 'Hide what changes' : 'What changes'}
                        title="What changes"
                        onClick={() => setOpen((value) => !value)}
                    >
                        ?
                    </button>
                )}
            </div>
            {open && details.length > 0 && (
                <ul className="llm-change-details" id={detailsId}>
                    {details.map((text) => <li key={text}>{text}</li>)}
                </ul>
            )}
            {note && <p className="llm-step-headline is-warn">{note}</p>}
            {said && <p className={`llm-step-headline${state === 'stale' || state === 'not_applied' ? ' is-fail' : ''}`}>{said}</p>}
            {undo && undo.state !== 'running' && undo.state !== 'done' && (
                <p className="llm-step-headline is-fail">{`Not undone: ${undo.reason}.`}</p>
            )}
            {undo && UNDO_TEXT[undo.state] && <p className="llm-step-headline">{UNDO_TEXT[undo.state]}</p>}
            {waiting && (
                <div className="llm-change-actions">
                    <button type="button" className="llm-primary llm-change-apply" onClick={() => onDecide?.(true)}>Apply</button>
                    <button type="button" className="llm-change-decline" onClick={() => onDecide?.(false)}>Decline</button>
                </div>
            )}
            {!waiting && (applied || pageLabel) && (
                <div className="llm-change-actions">
                    {applied && record?.id && (!undo || undo.state === 'failed') && (
                        <button
                            type="button"
                            className="ui-pill llm-change-undo"
                            onClick={() => onUndo?.(step)}
                            title="Put back the inputs the page had before this change"
                        >
                            Undo
                        </button>
                    )}
                    {pageLabel && onOpenPage && step.status !== 'running' && (
                        <button type="button" className="ui-pill" onClick={() => onOpenPage(step.page, null)} title={`Go to the ${pageLabel} page`}>
                            {`Open ${pageLabel}`}
                        </button>
                    )}
                </div>
            )}
        </li>
    );
};

export default ChangeCard;
