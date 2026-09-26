// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

import React from 'react';
import cx from './cx';

/**
 * Segmented toggle container.
 *
 * @param {string} [as='div'] - 'nav' for the app tabs.
 * @param {'frame'|'overlay'|'nav'} variant - card-header frame toggle, 3D
 *        viewport overlay pills, or the app navigation tabs.
 */
export const Segmented = ({ as = 'div', variant = 'frame', className, ...rest }) => {
    const Tag = as;
    return (
        <Tag className={cx('ui-seg', `ui-seg--${variant}`, className)} {...rest} />
    );
};

/**
 * One segment. Adds no `type`: pass it where the call site needs it.
 *
 * @param {boolean} [active]
 * @param {boolean} [overlay] - overlay-pill look (ui-seg__btn).
 * @param {boolean} [warm]    - warm active tint (overlay only).
 */
export const SegmentedButton = ({ active, overlay, warm, className, ...rest }) => (
    <button
        className={cx(overlay && 'ui-seg__btn', overlay && warm && 'ui-seg__btn--warm', active && 'is-active', className)}
        {...rest}
    />
);
