// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang
/* @vitest-environment jsdom */

// The shared Problems section: every error and warning on a page is listed at
// its top, nothing is dropped silently. Conditions hold while their message is
// set; one-off failures stay until dismissed, the same action succeeds, or
// another run opens; an error nothing caught lands on the page on screen.

import { afterEach, beforeEach, describe, expect, it, vi } from 'vitest';
import React, { act, useState } from 'react';
import { createRoot } from 'react-dom/client';
import { issueReducer, pageIssues, useIssue } from '../issueStore';
import { IssueList, IssueScope, IssueStoreProvider, IssueWatcher, PageIssues } from '../Issues';
import SaveMenu from '../SaveMenu';

const flush = () => act(async () => {
    for (let i = 0; i < 5; i += 1) await Promise.resolve();
});

describe('issueReducer', () => {
    const set = (state, key, message, extra = {}) => issueReducer(state, { type: 'set', page: 'p', key, message, severity: 'error', source: 'S', ...extra });
    const report = (state, message, extra = {}) => issueReducer(state, { type: 'report', page: 'p', message, severity: 'error', source: 'Save', ...extra });

    it('lists a condition while its message is set and removes it when cleared', () => {
        let state = set([], 'sites', 'Could not read the structure file.');
        expect(state).toHaveLength(1);
        expect(set(state, 'sites', 'Could not read the structure file.')).toBe(state);
        state = set(state, 'sites', null);
        expect(state).toEqual([]);
    });

    it('counts a repeated one-off failure instead of adding rows, keeping its row id', () => {
        let state = report([], 'Could not encode the figure');
        const { id } = state[0];
        state = report(state, 'Could not encode the figure');
        expect(state).toHaveLength(1);
        expect(state[0].count).toBe(2);
        expect(state[0].id).toBe(id);
    });

    it('keeps two reporters apart by key, and resolve clears only its own', () => {
        let state = report([], 'Could not encode the figure', { key: 'save:a' });
        state = report(state, 'Could not encode the figure', { key: 'save:b' });
        expect(state).toHaveLength(2);
        state = issueReducer(state, { type: 'resolve', page: 'p', key: 'save:b' });
        expect(state.map((issue) => issue.key)).toEqual(['save:a']);
    });

    it('dismissing a shown row dismisses every report it stands for', () => {
        let state = set([], 'structure-atoms-skipped', 'skipped 3 atoms', { severity: 'warning', source: 'Atoms skipped' });
        state = set(state, 'sites-atoms-skipped', 'skipped 3 atoms', { severity: 'warning', source: 'Atoms skipped' });
        const [shown] = pageIssues(state, 'p');
        state = issueReducer(state, { type: 'dismiss', id: shown.id });
        expect(pageIssues(state, 'p')).toEqual([]);
    });

    it('dismiss removes a one-off failure and hides a condition until its message changes', () => {
        let state = report(set([], 'sites', 'A'), 'B');
        const condition = state.find((issue) => issue.key === 'sites');
        const failure = state.find((issue) => issue.transient);
        state = issueReducer(state, { type: 'dismiss', id: failure.id });
        state = issueReducer(state, { type: 'dismiss', id: condition.id });
        expect(state).toHaveLength(1);
        expect(pageIssues(state, 'p')).toEqual([]);
        expect(pageIssues(set(state, 'sites', 'A'), 'p')).toEqual([]);
        expect(pageIssues(set(state, 'sites', 'A2'), 'p')).toHaveLength(1);
    });

    it('resolve and clearTransient drop one-off failures only; clearTransient spares `keep` pages', () => {
        let state = report(set([], 'sites', 'A'), 'B', { key: 'save:x' });
        expect(issueReducer(state, { type: 'resolve', page: 'p', key: 'save:x' })).toHaveLength(1);
        expect(issueReducer(state, { type: 'resolve', page: 'p', key: 'save:y' })).toBe(state);
        state = issueReducer(state, { type: 'report', page: 'autostog', message: 'bad upload', severity: 'error', source: null });
        state = issueReducer(state, { type: 'clearTransient', keep: ['autostog'] });
        expect(state.map((issue) => issue.key ?? issue.page)).toEqual(['sites', 'autostog']);
    });

    it('sync replaces a namespace and keeps unchanged entries (and their dismissal)', () => {
        let state = issueReducer([], { type: 'sync', page: 'p', namespace: 'files', items: [{ key: 'a', message: 'bad a' }, { key: 'b', message: 'bad b' }] });
        const [a] = state;
        state = issueReducer(state, { type: 'dismiss', id: a.id });
        const again = issueReducer(state, { type: 'sync', page: 'p', namespace: 'files', items: [{ key: 'a', message: 'bad a' }, { key: 'b', message: 'bad b' }] });
        expect(again).toBe(state);
        const fewer = issueReducer(state, { type: 'sync', page: 'p', namespace: 'files', items: [{ key: 'b', message: 'bad b' }, { key: 'c', message: null }] });
        expect(fewer.map((issue) => issue.message)).toEqual(['bad b']);
    });

    it('shows errors first, newest first, one row per distinct report, per page', () => {
        let state = set([], 'w', 'skipped 3 atoms', { severity: 'warning' });
        state = set(state, 'e1', 'first error');
        state = set(state, 'e2', 'second error');
        state = set(state, 'dup', 'skipped 3 atoms', { severity: 'warning' });
        state = issueReducer(state, { type: 'set', page: 'other', key: 'x', message: 'elsewhere', severity: 'error' });
        expect(pageIssues(state, 'p').map((issue) => issue.message)).toEqual(['second error', 'first error', 'skipped 3 atoms']);
        // The same message from two different figures is two problems.
        state = report(state, 'Could not encode the figure', { source: 'Save · F(Q)' });
        state = report(state, 'Could not encode the figure', { source: 'Save · G(r)' });
        expect(pageIssues(state, 'p').filter((issue) => issue.transient)).toHaveLength(2);
    });
});

describe('the section in a page', () => {
    let container;
    let root;
    beforeEach(() => {
        globalThis.IS_REACT_ACT_ENVIRONMENT = true;
        container = document.createElement('div');
        document.body.appendChild(container);
        root = createRoot(container);
    });
    afterEach(() => {
        act(() => root.unmount());
        container.remove();
        vi.restoreAllMocks();
    });

    const rows = () => [...container.querySelectorAll('.ui-issue')].map((row) => row.textContent);
    const inPage = (children, page = 'p', resetKey = 'run-1') => (
        <IssueStoreProvider>
            <IssueWatcher page={page} resetKey={resetKey} />
            <IssueScope page={page}>
                <PageIssues />
                {children}
            </IssueScope>
        </IssueStoreProvider>
    );

    const Condition = ({ message }) => {
        useIssue('sites', message, { source: 'Sites' });
        return null;
    };

    it('renders nothing visible without a problem, so the page layout is unchanged', () => {
        act(() => root.render(inPage(<Condition message={null} />)));
        expect(container.querySelector('.ui-issues')).toBeNull();
        // Only the two hidden live regions, empty, which take no space.
        const live = [...container.querySelectorAll('.ui-visually-hidden')];
        expect(live.map((node) => node.getAttribute('role'))).toEqual(['alert', 'status']);
        expect(live.every((node) => node.textContent === '')).toBe(true);
    });

    it('one save button succeeding does not clear another same-named button\'s failure', async () => {
        const Pair = () => (
            <>
                <div className="first"><SaveMenu onSave={() => Promise.reject(new Error('Could not encode the figure'))} name="G(r)" /></div>
                <div className="second"><SaveMenu onSave={() => Promise.resolve()} name="G(r)" /></div>
            </>
        );
        act(() => root.render(inPage(<Pair />)));
        await act(async () => container.querySelector('.first .ui-save__trigger').click());
        await flush();
        await act(async () => container.querySelector('.second .ui-save__trigger').click());
        await flush();
        expect(rows()).toEqual([expect.stringContaining('Could not encode the figure')]);
    });

    it('lists a condition with its source while it holds, and drops it when the component unmounts', () => {
        act(() => root.render(inPage(<Condition message="Request failed" />)));
        expect(rows()).toEqual([expect.stringContaining('Request failed')]);
        expect(container.querySelector('.ui-issue__source').textContent).toBe('Sites');
        expect(container.querySelector('.ui-issues').getAttribute('aria-label')).toBe('Problems on this page: 1 error');
        // The rows stay a list; the newest error is announced by a live region that is always there.
        expect(container.querySelector('.ui-issues__list').getAttribute('role')).toBeNull();
        expect(container.querySelector('[role="alert"]').textContent).toBe('Sites: Request failed');
        act(() => root.render(inPage(null)));
        expect(container.querySelector('.ui-issues')).toBeNull();
    });

    it('lists a failed save until the same save succeeds; outcomes that throw synchronously too', async () => {
        const outcomes = [Promise.reject(new Error('Could not encode the figure')), Promise.resolve()];
        const onSave = vi.fn(() => outcomes.shift());
        act(() => root.render(inPage(<SaveMenu onSave={onSave} name="3D view" />)));
        await act(async () => container.querySelector('.ui-save__trigger').click());
        await flush();
        expect(rows()).toEqual([expect.stringContaining('Could not encode the figure')]);
        expect(rows()[0]).toContain('Save · 3D view');
        await act(async () => container.querySelector('.ui-save__trigger').click());
        await flush();
        expect(rows()).toEqual([]);

        const throwing = () => { throw new Error('No canvas'); };
        act(() => root.render(inPage(<SaveMenu onSave={throwing} />)));
        await act(async () => container.querySelector('.ui-save__trigger').click());
        expect(rows()).toEqual([expect.stringContaining('No canvas')]);
    });

    it('outside a page scope, a failed save is logged, not dropped', async () => {
        const log = vi.spyOn(console, 'error').mockImplementation(() => {});
        act(() => root.render(<SaveMenu onSave={() => Promise.reject(new Error('nope'))} />));
        await act(async () => container.querySelector('.ui-save__trigger').click());
        await flush();
        expect(log).toHaveBeenCalledWith('Save failed:', expect.objectContaining({ message: 'nope' }));
    });

    it('lists an error nothing caught on the page on screen, and ignores browser noise', async () => {
        act(() => root.render(inPage(null, 'geometry')));
        await act(async () => {
            window.dispatchEvent(Object.assign(new Event('unhandledrejection'), { reason: new Error('worker exploded') }));
            window.dispatchEvent(Object.assign(new Event('error'), { message: 'ResizeObserver loop completed with undelivered notifications.' }));
        });
        expect(rows()).toEqual([expect.stringContaining('worker exploded')]);
        expect(rows()[0]).toContain('Unexpected error');
    });

    it('clears one-off failures when another run opens, and keeps conditions', async () => {
        const Page = () => {
            const [, force] = useState(0);
            return (
                <>
                    <Condition message="still broken" />
                    <SaveMenu onSave={() => Promise.reject(new Error('save failed'))} />
                    <button type="button" className="rerender" onClick={() => force((n) => n + 1)} />
                </>
            );
        };
        act(() => root.render(inPage(<Page />, 'p', 'run-1')));
        await act(async () => container.querySelector('.ui-save__trigger').click());
        await flush();
        expect(rows()).toHaveLength(2);
        act(() => root.render(inPage(<Page />, 'p', 'run-2')));
        expect(rows()).toEqual([expect.stringContaining('still broken')]);
    });

    it('is at most two lines collapsed: two problems, or the first and "N more"; a dismissed row goes', () => {
        const issues = [1, 2, 3, 4].map((n) => ({ id: n, severity: n < 3 ? 'error' : 'warning', source: null, message: `problem ${n}`, transient: true, count: n === 1 ? 3 : 1 }));
        act(() => root.render(<IssueList issues={issues.slice(0, 2)} />));
        expect(rows()).toHaveLength(2);
        expect(container.querySelector('.ui-issues__more')).toBeNull();

        const onDismiss = vi.fn();
        act(() => root.render(<IssueList issues={issues} onDismiss={onDismiss} />));
        expect(rows()).toHaveLength(1);
        expect(container.querySelector('.ui-issue__count').textContent).toBe('×3');
        const more = container.querySelector('.ui-issues__more');
        expect(more.textContent).toBe('3 more');
        act(() => more.click());
        expect(rows()).toHaveLength(4);
        expect(container.querySelector('.ui-issues__more').textContent).toBe('Show fewer');
        act(() => container.querySelector('.ui-issue__close').click());
        expect(onDismiss).toHaveBeenCalledWith(issues[0]);
        expect(container.querySelector('.ui-issues').getAttribute('aria-label')).toBe('Problems on this page: 2 errors, 2 warnings');
    });

    it('a page with warnings only is announced politely, not as an alert', () => {
        act(() => root.render(<IssueList issues={[{ id: 1, severity: 'warning', source: 'Atoms skipped', message: 'skipped', transient: false, count: 1 }]} />));
        expect(container.querySelector('.ui-issues').classList.contains('is-error')).toBe(false);
        expect(container.querySelector('[role="alert"]').textContent).toBe('');
        expect(container.querySelector('[role="status"]').textContent).toBe('Atoms skipped: skipped');
    });
});
