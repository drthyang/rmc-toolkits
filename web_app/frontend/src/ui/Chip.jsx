// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

import React from 'react';
import cx from './cx';

/**
 * Small pill-shaped readout.
 *
 * @param {'success'|'warn'|'danger'} [tone]
 * @param {boolean} [strong]   - heavier weight.
 * @param {boolean} [center]   - self-centered in a bottom-aligned row.
 * @param {boolean} [truncate] - ellipsize long values.
 */
const Chip = ({ tone, strong, center, truncate, className, ...rest }) => (
    <span
        className={cx(
            'ui-chip',
            strong && 'ui-chip--strong',
            center && 'ui-chip--center',
            truncate && 'ui-chip--truncate',
            tone && `ui-chip--${tone}`,
            className
        )}
        {...rest}
    />
);

/**
 * Element chip: a pill with a dot in the element color. Join chips with
 * BondDash inside a `ui-element-chain` span.
 *
 * @param {string}  color     - the element color (sets --chip).
 * @param {boolean} [central] - ring it (the central atom of a triplet).
 */
export const ElementChip = ({ color, central, className, style, children, ...rest }) => (
    <span
        className={cx('ui-element-chip', central && 'ui-element-chip--central', className)}
        style={color ? { '--chip': color, ...style } : style}
        {...rest}
    >
        <i className="ui-element-chip__dot" aria-hidden="true" />
        {children}
    </span>
);

/**
 * Bond dash: a short bar in a bond-role color. Its text (visually hidden,
 * default "–") keeps the chain readable as text: "Se–Ta–Se".
 *
 * @param {string}  color  - the bond color (sets --bond).
 * @param {boolean} [lead] - spacing after it, when it leads a chip's text.
 * @param {string}  [text='–'] - the hidden text.
 */
export const BondDash = ({ color, lead, text = '–', className, style, ...rest }) => (
    <span
        className={cx('ui-bond-dash', lead && 'ui-bond-dash--lead', className)}
        style={color ? { '--bond': color, ...style } : style}
        {...rest}
    >
        <span className="ui-visually-hidden">{text}</span>
    </span>
);

export default Chip;
