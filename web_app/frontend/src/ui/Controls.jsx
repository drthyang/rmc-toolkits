// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

import React from 'react';
import cx from './cx';

/**
 * Horizontal controls bar above a page's cards.
 *
 * @param {'default'|'dense'|'stacked'} [variant='default']
 * @param {boolean} [sub]    - a second bar tucked under the first.
 * @param {boolean} [footer] - a bar that closes the page.
 */
export const ControlsBar = ({ variant = 'default', sub, footer, className, ...rest }) => (
    <div
        className={cx(
            'ui-controls',
            variant !== 'default' && `ui-controls--${variant}`,
            sub && 'ui-controls--sub',
            footer && 'ui-controls--footer',
            className
        )}
        {...rest}
    />
);

/** A cluster of related controls that wraps as a unit (role="group"). */
export const ControlGroup = ({ label, className, ...rest }) => (
    <div className={cx('ui-control-group', className)} role="group" aria-label={label} {...rest} />
);

/**
 * One labeled control row: micro-label, the widget(s), an optional value.
 *
 * @param {string} [as='label'] - 'div' when the row holds more than one
 *                                interactive element (e.g. info badge + switch).
 * @param {React.ReactNode} label - micro-label content (text, info badge).
 * @param {React.ReactNode} [value] - readout after the widget.
 * @param {boolean} [valueWide] - wider readout column.
 */
export const Control = ({ as = 'label', label, value, valueWide, className, children, ...rest }) => {
    const Tag = as;
    return (
        <Tag className={cx('ui-control', className)} {...rest}>
            <span className="ui-control-label">{label}</span>
            {children}
            {value !== undefined && value !== null && (
                <span className={cx('ui-control-value', valueWide && 'ui-control-value--wide')}>{value}</span>
            )}
        </Tag>
    );
};

/**
 * Solid pill switch for a boolean option. Without `label` (and with `bare`)
 * it renders the pill alone, for rows whose label sits outside the <label>.
 */
export const Switch = ({ label, bare, checked, onChange, inputProps, className, ...rest }) => (
    <label className={cx('ui-control', 'ui-switch', bare && 'ui-switch--bare', className)} {...rest}>
        {label !== undefined && <span className="ui-control-label">{label}</span>}
        <input type="checkbox" checked={checked} onChange={onChange} {...inputProps} />
        <i className="ui-switch__track" aria-hidden="true" />
    </label>
);
