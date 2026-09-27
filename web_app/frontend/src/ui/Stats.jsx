// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

import React from 'react';
import cx from './cx';

/**
 * Stat rail: a titled cell on the left, labeled stat columns (<Stat>) beside.
 *
 * @param {React.ReactNode} heading - title cell content (title, source, badge).
 * @param {object} [headingProps]   - extra attributes for the <h2> (e.g. title).
 */
export const StatRail = ({ heading, headingProps, className, children, ...rest }) => (
    <section className={cx('ui-card', 'ui-stat-rail', className)} {...rest}>
        <h2 className="ui-stat-rail__title" {...headingProps}>{heading}</h2>
        <dl className="ui-stat-rail__stats">{children}</dl>
    </section>
);

/**
 * One stat column: <dt>{label}</dt><dd>{children}</dd>.
 *
 * @param {boolean} [end] - trailing group, flushed right.
 * @param {object} [dtProps] / [ddProps] - extra attributes for dt / dd.
 */
export const Stat = ({ label, end, dtProps, ddProps, className, children, ...rest }) => (
    <div className={cx('ui-stat', end && 'ui-stat--end', className)} {...rest}>
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
