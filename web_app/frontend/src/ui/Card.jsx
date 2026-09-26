// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

import React from 'react';
import cx from './cx';

/**
 * Card surface.
 *
 * @param {string}  [as='div'] - element (div / article / section / details).
 * @param {boolean} [clip]      - clip the content to the rounded corners.
 * @param {boolean} [lift]      - hover lift (stronger border + shadow).
 * @param {boolean} [roundEnds] - keep overflow visible (popovers escape); the
 *                                first / last children round the corners.
 * @param {'plot'|'bar'} [pad]  - inner padding preset.
 * @param {boolean} [note]      - small muted text card.
 */
export const Card = ({ as = 'div', clip, lift, roundEnds, pad, note, className, ...rest }) => {
    const Tag = as;
    return (
        <Tag
            className={cx(
                'ui-card',
                clip && 'ui-card--clip',
                roundEnds && 'ui-card--round-ends',
                lift && 'ui-card--lift',
                pad && `ui-card--pad-${pad}`,
                note && 'ui-card--note',
                className
            )}
            {...rest}
        />
    );
};

/**
 * Bar header of a card: title (+ help badge) on the left, then either a
 * `meta` node or an `actions` cluster on the right. Pass `children` instead of
 * the slots when the markup differs.
 *
 * @param {string}  [as='h3'] - 'div' where the header must not be a heading.
 * @param {boolean} [wrap]    - let the actions wrap under the title.
 * @param {boolean} [fixed]   - pinned height.
 */
export const CardHeader = ({ as = 'h3', wrap, fixed, title, help, meta, actions, children, className, ...rest }) => {
    const Tag = as;
    return (
        <Tag
            className={cx('ui-card__header', wrap && 'ui-card__header--wrap', fixed && 'ui-card__header--fixed', className)}
            {...rest}
        >
            {children !== undefined ? children : (
                <>
                    <span className="ui-card__label">{title}{help}</span>
                    {meta}
                    {actions !== undefined && <span className="ui-card__actions">{actions}</span>}
                </>
            )}
        </Tag>
    );
};

export const CardTitle = ({ as = 'h3', className, ...rest }) => {
    const Tag = as;
    return (
        <Tag className={cx('ui-card__title', className)} {...rest} />
    );
};

export const CardActions = ({ className, ...rest }) => (
    <span className={cx('ui-card__actions', className)} {...rest} />
);

export const CardMeta = ({ fixed, className, ...rest }) => (
    <span className={cx('ui-card__meta', fixed && 'ui-card__meta--fixed', className)} {...rest} />
);

export const CardNote = ({ emph, className, ...rest }) => (
    <div className={cx('ui-card__note', emph && 'ui-card__note--emph', className)} {...rest} />
);
