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

export default Chip;
