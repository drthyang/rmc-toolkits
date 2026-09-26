// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

import React from 'react';
import cx from './cx';
import { IconButton } from './Buttons';

/**
 * Message banner.
 *
 * @param {string} [as='div'] - 'p' for a single paragraph.
 * @param {'danger'|'neutral'|'caution'|'danger-light'} tone
 * @param {boolean} [sm] / [gapLg] / [flush] / [inline] - size and spacing presets.
 * @param {Function} [onDismiss] - adds a close (×) button beside the message.
 */
export const Banner = ({ as = 'div', tone, sm, gapLg, flush, inline, onDismiss, className, children, ...rest }) => {
    const Tag = as;
    return (
        <Tag
            className={cx(
                'ui-banner',
                tone && `ui-banner--${tone}`,
                inline && 'ui-banner--inline',
                sm && 'ui-banner--sm',
                gapLg && 'ui-banner--gap-lg',
                flush && 'ui-banner--flush',
                onDismiss && 'ui-banner--dismissible',
                className
            )}
            {...rest}
        >
            {onDismiss ? (
                <>
                    <span>{children}</span>
                    <IconButton variant="close" onClick={onDismiss} aria-label="Close notification" title="Close">
                        &times;
                    </IconButton>
                </>
            ) : children}
        </Tag>
    );
};

/** Dashed hint box (what to do next). */
export const Hint = ({ className, ...rest }) => (
    <p className={cx('ui-hint', className)} {...rest} />
);

/** Empty-state card; `fill` grows it to the free height and centers the text. */
export const EmptyState = ({ fill, className, ...rest }) => (
    <div className={cx('ui-empty', fill && 'ui-empty--fill', className)} {...rest} />
);
