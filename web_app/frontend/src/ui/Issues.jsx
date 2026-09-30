// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

import React, { useContext, useEffect, useReducer, useRef, useState } from 'react';
import cx from './cx';
import { IconButton } from './Buttons';
import {
    IssueDispatchContext,
    IssuePageContext,
    IssueStateContext,
    issueReducer,
    usePageIssues
} from './issueStore';

// Collapsed, the section is at most two lines tall: two problems show both;
// more show the first (the newest error) and "N more".
const COLLAPSED_LINES = 2;

/** The app's one problems store (see issueStore.js). */
export const IssueStoreProvider = ({ children }) => {
    const [issues, dispatch] = useReducer(issueReducer, []);
    return (
        <IssueDispatchContext.Provider value={dispatch}>
            <IssueStateContext.Provider value={issues}>
                {children}
            </IssueStateContext.Provider>
        </IssueDispatchContext.Provider>
    );
};

/** Everything inside reports to page `page`'s Problems section. */
export const IssueScope = ({ page, children }) => (
    <IssuePageContext.Provider value={page}>{children}</IssuePageContext.Provider>
);

/**
 * The Problems section: one compact row per problem (severity mark, source,
 * message, a ×N count for a repeated failure, a close button). At most two
 * lines show: two problems, or the first and "N more". Renders nothing
 * visible when there is no problem (only its two hidden live regions).
 */
export const IssueList = ({ issues, onDismiss, className }) => {
    const [expanded, setExpanded] = useState(false);
    const errors = issues.filter((issue) => issue.severity === 'error');
    const warnings = issues.filter((issue) => issue.severity !== 'error');
    // Announcements go through two always-mounted, visually hidden live regions
    // whose roles never change (the list keeps its list semantics): the newest
    // error (assertive) and the newest warning (polite).
    const announce = (issue) => (issue ? `${issue.source ? `${issue.source}: ` : ''}${issue.message}` : '');
    const live = (
        <>
            <span className="ui-visually-hidden" role="alert">{announce(errors[0])}</span>
            <span className="ui-visually-hidden" role="status">{announce(warnings[0])}</span>
        </>
    );
    if (!issues.length) return live;
    const shown = expanded || issues.length <= COLLAPSED_LINES ? issues : issues.slice(0, COLLAPSED_LINES - 1);
    const hidden = issues.length - shown.length;
    const label = [
        errors.length && `${errors.length} error${errors.length === 1 ? '' : 's'}`,
        warnings.length && `${warnings.length} warning${warnings.length === 1 ? '' : 's'}`
    ].filter(Boolean).join(', ');
    return (
        <>
            {live}
            <section
                className={cx('ui-issues', errors.length > 0 && 'is-error', className)}
                aria-label={`Problems on this page: ${label}`}
            >
                <ul className="ui-issues__list">
                    {shown.map((issue) => (
                        <li key={issue.id} className={cx('ui-issue', `ui-issue--${issue.severity}`)}>
                            <span className="ui-issue__mark" aria-hidden="true" />
                            <span className="ui-visually-hidden">{issue.severity === 'error' ? 'Error' : 'Warning'}: </span>
                            {issue.source && <span className="ui-issue__source">{issue.source}</span>}
                            <span className="ui-issue__text">{issue.message}</span>
                            {issue.count > 1 && (
                                <span className="ui-issue__count" title={`Happened ${issue.count} times`}>×{issue.count}</span>
                            )}
                            {onDismiss && (
                                <IconButton
                                    variant="close"
                                    className="ui-issue__close"
                                    onClick={() => onDismiss(issue)}
                                    aria-label="Dismiss this problem"
                                    title={issue.transient ? 'Dismiss' : 'Hide until it changes'}
                                >
                                    &times;
                                </IconButton>
                            )}
                        </li>
                    ))}
                </ul>
                {(hidden > 0 || expanded) && issues.length > COLLAPSED_LINES && (
                    <button
                        type="button"
                        className="ui-issues__more"
                        onClick={() => setExpanded((value) => !value)}
                        aria-expanded={expanded}
                    >
                        {expanded ? 'Show fewer' : `${hidden} more`}
                    </button>
                )}
            </section>
        </>
    );
};

/** This page's Problems section, where the page places it (its first child). */
export const PageIssues = ({ className }) => {
    const issues = usePageIssues();
    const dispatch = useContext(IssueDispatchContext);
    return (
        <IssueList
            issues={issues}
            className={className}
            onDismiss={dispatch ? (issue) => dispatch({ type: 'dismiss', id: issue.id }) : undefined}
        />
    );
};

// Browser noise that is not a failure of the app.
const IGNORED = [/ResizeObserver loop/i];
const isIgnored = (reason) => {
    if (!reason) return true;
    if (reason.name === 'AbortError' || reason.name === 'CanceledError') return true;
    const message = reason.message ?? String(reason);
    return IGNORED.some((pattern) => pattern.test(message));
};

/**
 * App-level watcher: an error nothing caught (an exception in an event
 * handler, a rejected promise) is listed in the active page's Problems
 * section instead of vanishing into the console; and when `resetKey` (the
 * open run) changes, one-off failures of the previous run are cleared, except
 * on the pages in `keep` (those that do not depend on the run).
 */
const NO_PAGES = [];
export const IssueWatcher = ({ page, resetKey, keep = NO_PAGES }) => {
    const dispatch = useContext(IssueDispatchContext);
    const pageRef = useRef(page);
    useEffect(() => {
        pageRef.current = page;
    }, [page]);

    useEffect(() => {
        if (!dispatch) return undefined;
        const report = (reason) => {
            if (isIgnored(reason)) return;
            const message = reason?.message || String(reason);
            dispatch({ type: 'report', page: pageRef.current, severity: 'error', source: 'Unexpected error', message });
        };
        const onError = (event) => report(event.error ?? (event.message ? { message: event.message } : null));
        const onRejection = (event) => report(event.reason);
        window.addEventListener('error', onError);
        window.addEventListener('unhandledrejection', onRejection);
        return () => {
            window.removeEventListener('error', onError);
            window.removeEventListener('unhandledrejection', onRejection);
        };
    }, [dispatch]);

    const keepRef = useRef(keep);
    useEffect(() => {
        keepRef.current = keep;
    }, [keep]);
    const firstRun = useRef(true);
    useEffect(() => {
        if (firstRun.current) {
            firstRun.current = false;
            return;
        }
        if (dispatch) dispatch({ type: 'clearTransient', keep: keepRef.current });
    }, [dispatch, resetKey]);

    return null;
};
