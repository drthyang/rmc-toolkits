// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

import React from 'react';
import cx from './cx';

/**
 * Pill button.
 *
 * @param {'sm'|'md'} [size='sm'] - 'md' is the taller pill with an active state.
 * @param {boolean} [tint]   - accent-tinted pill (e.g. Reset zoom).
 * @param {boolean} [active] - selected state (md).
 */
export const Pill = ({ size = 'sm', tint, active, className, ...rest }) => (
    <button
        type="button"
        className={cx(tint ? 'ui-pill-tint' : size === 'md' ? 'ui-pill-md' : 'ui-pill', active && 'is-active', className)}
        {...rest}
    />
);

/** Card-header tool chip (Reset view; `axes` for the a b c toggle). */
export const ToolButton = ({ axes, active, className, ...rest }) => (
    <button
        type="button"
        className={cx('ui-tool-btn', axes && 'ui-tool-btn--axes', active && 'is-active', className)}
        {...rest}
    />
);

/** Round × button: 'close' (notifications) or 'remove' (file chips). */
export const IconButton = ({ variant = 'close', className, ...rest }) => (
    <button type="button" className={cx('ui-icon-btn', `ui-icon-btn--${variant}`, className)} {...rest} />
);

/** Solid accent button. Adds no `type`: the caller passes it. */
export const PrimaryButton = ({ outlined, className, ...rest }) => (
    <button className={cx('ui-btn-primary', outlined && 'ui-btn-primary--outlined', className)} {...rest} />
);
