// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// The page problems store behind the shared "Problems" section (Issues.jsx).
//
// Every workspace page lists what went wrong on it in one compact section at
// its top: errors and warnings, both the ones that hold while a condition
// lasts (a site table that failed to load, atoms the parser skipped) and the
// ones an action raised once (a figure that could not be saved). Nothing is
// dropped silently. One store for the app (IssueStoreProvider); each page is
// a scope (IssueScope page="…") that its components report into.
//
//   useIssue(key, message, { severity, source })  a condition: listed while
//                                                  `message` is set, removed
//                                                  when it clears or unmounts
//   useIssueSet(namespace, items)                  a list of conditions (one per
//                                                  plot file, say) kept in sync
//   useReportIssue()                               a one-off failure: report()
//                                                  lists it until dismissed,
//                                                  another run opens, or the
//                                                  same action succeeds
//   useResolveIssue()                              resolve(key): the action that
//                                                  reported under `key` succeeded
//
// Hooks read only the stable dispatch, so reporting never re-renders a page.

import { createContext, useCallback, useContext, useEffect } from 'react';

export const IssueDispatchContext = createContext(null);
export const IssueStateContext = createContext([]);
export const IssuePageContext = createContext(null);

let nextId = 0;

// A message as display text, or null for nothing to show.
const text = (message) => {
    if (message == null || message === false) return null;
    const value = typeof message === 'string' ? message : (message?.message ?? String(message));
    return value.trim() ? value : null;
};

const sameContent = (issue, action) => issue.message === action.message
    && issue.severity === action.severity
    && issue.source === action.source;

export const issueReducer = (state, action) => {
    switch (action.type) {
        case 'set': {
            // A condition: one entry per (page, key); cleared with a null message.
            const index = state.findIndex((issue) => issue.page === action.page && issue.key === action.key);
            if (!action.message) {
                return index < 0 ? state : state.filter((_, i) => i !== index);
            }
            if (index >= 0 && sameContent(state[index], action)) return state;
            nextId += 1;
            const entry = {
                id: nextId,
                page: action.page,
                key: action.key,
                severity: action.severity,
                source: action.source,
                message: action.message,
                transient: false,
                dismissed: false,
                count: 1
            };
            return index < 0 ? [...state, entry] : state.map((issue, i) => (i === index ? entry : issue));
        }
        case 'sync': {
            // A namespace's conditions replaced by a list; unchanged entries keep
            // their id (and a dismissal), so a re-render does not re-show them.
            const prefix = `${action.namespace}\u0000`;
            const inNamespace = (issue) => issue.page === action.page && !issue.transient
                && typeof issue.key === 'string' && issue.key.startsWith(prefix);
            const current = state.filter(inNamespace);
            const next = action.items.filter((item) => text(item.message)).map((item) => {
                const key = `${prefix}${item.key}`;
                const message = text(item.message);
                const severity = item.severity ?? 'error';
                const source = item.source ?? null;
                const existing = current.find((issue) => issue.key === key);
                if (existing && sameContent(existing, { message, severity, source })) return existing;
                nextId += 1;
                return { id: nextId, page: action.page, key, severity, source, message, transient: false, dismissed: false, count: 1 };
            });
            if (next.length === current.length && next.every((issue, i) => issue === current[i])) return state;
            return [...state.filter((issue) => !inNamespace(issue)), ...next];
        }
        case 'report': {
            // A one-off failure: the same failure again (same reporter `key` when
            // given, e.g. one save button; else the same content) counts up
            // instead of adding a row, keeping its id so the row is not remounted
            // (an error thrown every frame does not churn the list), and comes
            // back if it had been dismissed.
            const same = (issue) => issue.transient && issue.page === action.page && sameContent(issue, action)
                && (action.key == null || issue.key === action.key);
            const index = state.findIndex(same);
            if (index >= 0) {
                const bumped = { ...state[index], count: state[index].count + 1, dismissed: false };
                return state.map((issue, i) => (i === index ? bumped : issue));
            }
            nextId += 1;
            return [...state, {
                id: nextId,
                page: action.page,
                key: action.key ?? null,
                severity: action.severity,
                source: action.source,
                message: action.message,
                transient: true,
                dismissed: false,
                count: 1
            }];
        }
        case 'resolve': {
            // The action succeeded: its earlier one-off failures no longer apply.
            // Matched by the reporter's key (one save button), not by its label,
            // which two charts can share.
            const resolved = (issue) => issue.transient && issue.page === action.page && issue.key === action.key;
            return state.some(resolved) ? state.filter((issue) => !resolved(issue)) : state;
        }
        case 'dismiss': {
            // A shown row can stand for several reports of the same thing (see
            // pageIssues): dismiss them all. A one-off failure goes; a condition
            // is hidden until its message changes.
            const target = state.find((issue) => issue.id === action.id);
            if (!target) return state;
            const group = (issue) => issue.page === target.page && issue.severity === target.severity
                && (issue.source ?? '') === (target.source ?? '') && issue.message === target.message;
            return state
                .filter((issue) => !(group(issue) && issue.transient))
                .map((issue) => (group(issue) ? { ...issue, dismissed: true } : issue));
        }
        case 'clearTransient': {
            // Another run opened: last run's one-off failures no longer apply
            // (except on pages that do not depend on the run, `keep`).
            const keep = action.keep ?? [];
            const cleared = (issue) => issue.transient && !keep.includes(issue.page);
            return state.some(cleared) ? state.filter((issue) => !cleared(issue)) : state;
        }
        default:
            return state;
    }
};

/**
 * List a condition in this page's Problems section while `message` is set.
 * `key` is unique within the page; severity 'error' (default) or 'warning';
 * `source` names what it concerns ("Sites", "Save · F(Q)").
 */
export const useIssue = (key, message, { severity = 'error', source = null } = {}) => {
    const dispatch = useContext(IssueDispatchContext);
    const page = useContext(IssuePageContext);
    const value = text(message);
    useEffect(() => {
        if (!dispatch || !page) return undefined;
        dispatch({ type: 'set', page, key, message: value, severity, source });
        return () => dispatch({ type: 'set', page, key, message: null });
    }, [dispatch, page, key, value, severity, source]);
};

/**
 * Keep a list of conditions — [{ key, message, severity?, source? }] — in sync
 * under `namespace`: entries with a message are listed, the rest removed.
 * Pass a memoised array.
 */
export const useIssueSet = (namespace, items) => {
    const dispatch = useContext(IssueDispatchContext);
    const page = useContext(IssuePageContext);
    useEffect(() => {
        if (!dispatch || !page) return undefined;
        dispatch({ type: 'sync', page, namespace, items: items || [] });
        return undefined;
    }, [dispatch, page, namespace, items]);
    useEffect(() => () => {
        if (dispatch && page) dispatch({ type: 'sync', page, namespace, items: [] });
    }, [dispatch, page, namespace]);
};

/**
 * A function that lists a one-off failure in this page's Problems section —
 * report({ message, severity?, source? }) — or null outside a page scope, so a
 * component can fall back to showing it itself.
 */
export const useReportIssue = () => {
    const dispatch = useContext(IssueDispatchContext);
    const page = useContext(IssuePageContext);
    const report = useCallback(({ message, severity = 'error', source = null, key = null }) => {
        const value = text(message);
        if (value) dispatch({ type: 'report', page, message: value, severity, source, key });
    }, [dispatch, page]);
    return dispatch && page ? report : null;
};

/** A function that clears this page's one-off failures reported under `key` (that action succeeded), or null. */
export const useResolveIssue = () => {
    const dispatch = useContext(IssueDispatchContext);
    const page = useContext(IssuePageContext);
    const resolve = useCallback((key) => dispatch({ type: 'resolve', page, key }), [dispatch, page]);
    return dispatch && page ? resolve : null;
};

/**
 * The shown problems of one page: errors first, newest first. Two reports of
 * the same thing (same severity, source and message: the skipped-atoms warning
 * of both the model card and the site table) show as one row.
 */
export const pageIssues = (issues, page) => {
    const seen = new Set();
    return issues
        .filter((issue) => issue.page === page && !issue.dismissed)
        .sort((a, b) => ((a.severity === 'error' ? 0 : 1) - (b.severity === 'error' ? 0 : 1)) || (b.id - a.id))
        .filter((issue) => {
            const signature = `${issue.severity}\u0000${issue.source ?? ''}\u0000${issue.message}`;
            if (seen.has(signature)) return false;
            seen.add(signature);
            return true;
        });
};

export const usePageIssues = () => {
    const issues = useContext(IssueStateContext);
    const page = useContext(IssuePageContext);
    return pageIssues(issues, page);
};
