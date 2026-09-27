// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

import React, { useLayoutEffect, useRef } from 'react';
import cx from './cx';
import { planStatRows } from './statRows';

// Stats that do not all fit on one line would wrap into ragged lines, each
// starting indented behind a stray divider. StatRail measures the natural
// width of the runs of stats between band starts (<Stat band>) in the plain
// flex row — where every stat has its natural width, so the answer does not
// depend on the layout it selects — and planStatRows picks the layout: the
// one-line row, `is-banded` (runs packed into rows, broken before a band
// start only where needed: data-break-a / data-break-b), or `is-wrapped` (an
// aligned grid, for a rail without bands). ui.css draws each. On a phone the
// list is always a grid; nothing is measured there.
const useWrapState = (listRef) => {
    useLayoutEffect(() => {
        const list = listRef.current;
        if (!list || typeof window === 'undefined') return undefined;
        const update = () => {
            list.classList.remove('is-wrapped', 'is-banded');
            list.removeAttribute('data-break-a');
            list.removeAttribute('data-break-b');
            const style = window.getComputedStyle(list);
            if (style.display !== 'flex' || list.children.length < 2) return;
            const available = list.clientWidth
                - (parseFloat(style.paddingLeft) || 0)
                - (parseFloat(style.paddingRight) || 0)
                + 0.5;
            const runs = [0];
            for (const item of list.children) {
                if (item.classList.contains('ui-stat--band') && item !== list.firstElementChild && runs.length < 3) {
                    runs.push(0);
                }
                runs[runs.length - 1] += item.getBoundingClientRect().width;
            }
            const plan = planStatRows(runs, available);
            if (plan.mode === 'wrapped') list.classList.add('is-wrapped');
            if (plan.mode === 'banded') {
                list.classList.add('is-banded');
                if (plan.breakA) list.setAttribute('data-break-a', '');
                if (plan.breakB) list.setAttribute('data-break-b', '');
            }
        };
        update();
        if (typeof ResizeObserver === 'undefined') return undefined;
        const observer = new ResizeObserver(update);
        observer.observe(list);
        return () => observer.disconnect();
    });
};

/**
 * Stat rail: a titled cell on the left, labeled stat columns (<Stat>) beside.
 * When the stats do not fit on one line they form an aligned grid instead.
 *
 * @param {React.ReactNode} heading - title cell content (title, source, badge).
 * @param {object} [headingProps]   - extra attributes for the <h2> (e.g. title).
 */
export const StatRail = ({ heading, headingProps, className, children, ...rest }) => {
    const listRef = useRef(null);
    useWrapState(listRef);
    return (
        <section className={cx('ui-card', 'ui-stat-rail', className)} {...rest}>
            <h2 className="ui-stat-rail__title" {...headingProps}>{heading}</h2>
            <dl className="ui-stat-rail__stats" ref={listRef}>{children}</dl>
        </section>
    );
};

/**
 * One stat column: <dt>{label}</dt><dd>{children}</dd>.
 *
 * @param {boolean} [end]  - trailing group, flushed right.
 * @param {boolean} [band] - starts a band: a new row when the rail's stats
 *                           do not fit on one line (at most two per rail).
 * @param {object} [dtProps] / [ddProps] - extra attributes for dt / dd.
 */
export const Stat = ({ label, end, band, dtProps, ddProps, className, children, ...rest }) => (
    <div className={cx('ui-stat', end && 'ui-stat--end', band && 'ui-stat--band', className)} {...rest}>
        <dt {...dtProps}>{label}</dt>
        <dd {...ddProps}>{children}</dd>
    </div>
);

/**
 * Readout tile with a status edge.
 *
 * @param {'good'|'warn'|'bad'} [tone]
 */
export const StatCard = ({ tone, label, value, sub, className, ...rest }) => (
    <div className={cx('ui-stat-card', tone && `is-${tone}`, className)} {...rest}>
        <span className="ui-stat-card__label">{label}</span>
        <span className="ui-stat-card__value">{value}</span>
        <span className="ui-stat-card__sub">{sub}</span>
    </div>
);

/**
 * KPI rail: headline results as a <dl> of tiles inside a card, under its
 * header. The caller adds role / aria-live / aria-label.
 */
export const KpiRail = ({ className, ...rest }) => (
    <dl className={cx('ui-kpis', className)} {...rest} />
);

/**
 * One KPI tile. Without a value it shows "—" and keeps its sub line (a
 * non-breaking space), so the rail keeps its height before the first result.
 *
 * @param {React.ReactNode} label
 * @param {React.ReactNode} [value] - null / undefined → "—".
 * @param {React.ReactNode} [unit]  - after the value, past a thin space.
 * @param {React.ReactNode} [sub]   - one muted line under the value.
 */
export const Kpi = ({ label, value, unit, sub, className, ...rest }) => {
    const empty = value === null || value === undefined;
    return (
        <div className={cx('ui-kpi', className)} {...rest}>
            <dt>{label}</dt>
            <dd>
                <span className={cx('ui-kpi__value', empty && 'is-empty')}>
                    {empty ? '—' : value}
                    {!empty && unit ? <span className="ui-kpi__unit">{'\u2009'}{unit}</span> : null}
                </span>
                <span className="ui-kpi__sub">{empty || sub === undefined || sub === null ? '\u00a0' : sub}</span>
            </dd>
        </div>
    );
};
